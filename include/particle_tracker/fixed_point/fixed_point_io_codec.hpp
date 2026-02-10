// fixed_point_io_codec.hpp
// PKD1 "record-stream" compressed helpers (zstd or gzip/zlib).
// Requires pk_config.hpp to be included first (or compile flags defined).

// NOTE ON SIZE LIMITS AND FORMAT CONSTRAINTS
//
// * gzip/zlib path:
//   - zlib APIs (deflate, crc32) use `uInt` for input sizes (typically 32-bit).
//   - A single gzip member therefore cannot safely process >4 GiB of
//     uncompressed data in one call.
//   - This implementation hard-fails on overflow instead of silently
//     truncating and corrupting the stream.
//   - Supporting >4 GiB would require explicit chunked deflate + rolling CRC.
//
// * zstd path:
//   - Zstd compression/decompression APIs use `size_t` and are not subject
//     to the 4 GiB zlib limitation.
//   - Zstd skippable frame headers store payload size as uint32_t by format;
//     this is safe here because FileHeaderPKD1 is small.
//   - Decompression accumulates the full output in memory; no hard size
//     limit is enforced for corrupted or untrusted inputs.

#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <cstring>

#include "pk_config.hpp"
#include "fixed_point_core.hpp"
#include "fixed_point_pkd1_io.hpp" // for write_all/read_all + FileHeaderPKD1

#ifdef PK_IO_USE_ZSTD
    #include <zstd.h>     // link with -lzstd
#endif
#ifdef PK_IO_USE_ZLIB
    #include <zlib.h>     // link with -lz
#endif

// -------------------- ZSTD implementation --------------------
#ifdef PK_IO_USE_ZSTD

inline void write_records_zstd(const std::string& path,
                               const std::vector<PackedParticle>& rec,
                               std::uint64_t released_count,
                               int level = 3,
                               bool append = false)
{
    if (!append)
    {
        std::ofstream os(path, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!os) throw std::runtime_error("write_records_zstd: cannot open " + path);

        FileHeaderPKD1 h{};
        h.reserved0 = released_count;
        h.reserved1 = static_cast<uint64_t>(rec.size());

        // skippable frame magic (0x184D2A50..57)
        const uint32_t magic = 0x184D2A50u;
        const uint32_t size  = static_cast<uint32_t>(sizeof(FileHeaderPKD1));
        os.write(reinterpret_cast<const char*>(&magic), 4);
        os.write(reinterpret_cast<const char*>(&size),  4);
        write_header(os, h);
        os.flush();
    }
    else
    {
        // bump reserved1 inside skippable frame at file start
        std::fstream hs(path, std::ios::binary | std::ios::in | std::ios::out);
        if (!hs) throw std::runtime_error("write_records_zstd: cannot reopen for header " + path);

        uint32_t magic = 0, size = 0;
        hs.read(reinterpret_cast<char*>(&magic), 4);
        hs.read(reinterpret_cast<char*>(&size),  4);
        if ((magic & 0xFFFFFFF0u) != 0x184D2A50u || size < sizeof(FileHeaderPKD1))
            throw std::runtime_error("write_records_zstd: no skippable header found");

        FileHeaderPKD1 h = read_header(hs);
        h.reserved1 += static_cast<uint64_t>(rec.size());

        hs.seekp(8, std::ios::beg); // after magic+size
        write_header(hs, h);
        hs.flush();
    }

    if (rec.empty()) return;

    std::ofstream os(path, std::ios::binary | std::ios::out | std::ios::app);
    if (!os) throw std::runtime_error("write_records_zstd: cannot open for append " + path);

    const size_t src_bytes = rec.size() * sizeof(PackedParticle);
    const size_t bound     = ZSTD_compressBound(src_bytes);
    std::vector<char> cbuf(bound);

    const size_t rc = ZSTD_compress(cbuf.data(), bound, rec.data(), src_bytes, level);
    if (ZSTD_isError(rc))
        throw std::runtime_error("write_records_zstd: " + std::string(ZSTD_getErrorName(rc)));

    write_all(os, cbuf.data(), rc);
}

inline std::vector<PackedParticle> read_records_zstd(const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("read_records_zstd: cannot open " + path);

    std::streampos start = is.tellg();
    uint32_t magic = 0, size = 0;
    is.read(reinterpret_cast<char*>(&magic), 4);
    is.read(reinterpret_cast<char*>(&size),  4);

    if (is && ((magic & 0xFFFFFFF0u) == 0x184D2A50u))
    {
        if (size < sizeof(FileHeaderPKD1))
            throw std::runtime_error("read_records_zstd: skippable too small");
        FileHeaderPKD1 h{};
        read_all(is, &h);
        if (h.magic[0] != 'P' || h.magic[1] != 'K' || h.magic[2] != 'D')
            throw std::runtime_error("read_records_zstd: bad PKD header");
        if (size > sizeof(FileHeaderPKD1))
            is.seekg(static_cast<std::streamoff>(size - sizeof(FileHeaderPKD1)), std::ios::cur);
    }
    else
    {
        // backward compat: old files had raw PKD1 header
        is.clear();
        is.seekg(start);
        (void)read_header(is);
    }

    ZSTD_DCtx* dctx = ZSTD_createDCtx();
    if (!dctx) throw std::runtime_error("read_records_zstd: createDCtx failed");

    const size_t in_chunk  = ZSTD_DStreamInSize();
    const size_t out_chunk = ZSTD_DStreamOutSize();
    std::vector<char> inbuf(in_chunk), outbuf(out_chunk), raw;

    ZSTD_inBuffer  in  { inbuf.data(), 0, 0 };
    ZSTD_outBuffer out { outbuf.data(), outbuf.size(), 0 };

    while (true)
    {
        if (in.pos == in.size)
        {
            is.read(inbuf.data(), static_cast<std::streamsize>(inbuf.size()));
            in.size = static_cast<size_t>(is.gcount());
            in.pos  = 0;
            if (in.size == 0) break;
        }

        out.pos = 0;
        size_t ret = ZSTD_decompressStream(dctx, &out, &in);
        if (ZSTD_isError(ret))
        {
            ZSTD_freeDCtx(dctx);
            throw std::runtime_error("read_records_zstd: " + std::string(ZSTD_getErrorName(ret)));
        }
        if (out.pos) raw.insert(raw.end(), outbuf.data(), outbuf.data() + out.pos);
    }

    ZSTD_freeDCtx(dctx);

    if (raw.size() % sizeof(PackedParticle) != 0)
        throw std::runtime_error("read_records_zstd: size mismatch");

    std::vector<PackedParticle> rec(raw.size() / sizeof(PackedParticle));
    std::memcpy(rec.data(), raw.data(), raw.size());
    return rec;
}

#endif // PK_IO_USE_ZSTD

// -------------------- ZLIB/GZIP implementation --------------------
#ifdef PK_IO_USE_ZLIB

namespace gzdetail
{
    inline uint16_t get_le16(std::istream& is)
    {
        unsigned char b[2]{};
        if (!is.read(reinterpret_cast<char*>(b), 2))
            throw std::runtime_error("get_le16: read failed");
        return static_cast<uint16_t>(b[0]) | (static_cast<uint16_t>(b[1]) << 8);
    }

    inline uint32_t get_le32(std::istream& is)
    {
        unsigned char b[4]{};
        if (!is.read(reinterpret_cast<char*>(b), 4))
            throw std::runtime_error("get_le32: read failed");
        return static_cast<uint32_t>(b[0])
             | (static_cast<uint32_t>(b[1]) << 8)
             | (static_cast<uint32_t>(b[2]) << 16)
             | (static_cast<uint32_t>(b[3]) << 24);
    }

    inline void put_le16(std::ostream& os, uint16_t v)
    {
        unsigned char b[2] = {
            static_cast<unsigned char>(v & 0xFFu),
            static_cast<unsigned char>((v >> 8) & 0xFFu)
        };
        os.write(reinterpret_cast<const char*>(b), 2);
    }

    inline void put_le32(std::ostream& os, uint32_t v)
    {
        unsigned char b[4] = {
            static_cast<unsigned char>(v & 0xFFu),
            static_cast<unsigned char>((v >> 8) & 0xFFu),
            static_cast<unsigned char>((v >> 16) & 0xFFu),
            static_cast<unsigned char>((v >> 24) & 0xFFu)
        };
        os.write(reinterpret_cast<const char*>(b), 4);
    }

    // gzip header with FEXTRA containing FileHeaderPKD1
    inline void write_gzip_header(std::ostream& os, const FileHeaderPKD1& h, bool with_extra)
    {
        const uint8_t ID1=0x1f, ID2=0x8b, CM=8;
        uint8_t FLG = with_extra ? 0x04 : 0x00; // FEXTRA

        os.put(ID1); os.put(ID2); os.put(CM); os.put(FLG);
        put_le32(os, 0);          // MTIME
        os.put(0);                // XFL
        unsigned char os_byte=0xFF;
        os.write(reinterpret_cast<const char*>(&os_byte), 1);

        if (with_extra)
        {
            const uint16_t sublen = static_cast<uint16_t>(sizeof(FileHeaderPKD1));
            const uint16_t xlen   = static_cast<uint16_t>(2 + 2 + sublen);
            put_le16(os, xlen);

            os.put('P'); os.put('K');      // subfield ID
            put_le16(os, sublen);
            write_header(os, h);
        }
    }

    inline void bump_gzip_extra_reserved1(const std::string& path, uint64_t add)
    {
        std::fstream hs(path, std::ios::binary | std::ios::in | std::ios::out);
        if (!hs) throw std::runtime_error("write_records_gzip: cannot reopen " + path);
    
        unsigned char hdr[10];
        hs.read(reinterpret_cast<char*>(hdr), 10);
        if (hs.gcount()!=10 || hdr[0]!=0x1f || hdr[1]!=0x8b || hdr[2]!=8)
        throw std::runtime_error("write_records_gzip: not a gzip file");
    
        const uint8_t flg = hdr[3];
        if ((flg & 0x04) == 0) return; // no FEXTRA
    
        const uint16_t xlen = gzdetail::get_le16(hs);
        if (xlen < 4)
        throw std::runtime_error("write_records_gzip: corrupt FEXTRA (xlen < 4)");
    
        char si1 = 0, si2 = 0;
        hs.get(si1); hs.get(si2);
        const uint16_t len = gzdetail::get_le16(hs);
    
        if (si1!='P' || si2!='K' || len < sizeof(FileHeaderPKD1))
        throw std::runtime_error("write_records_gzip: PKD extra not found");
    
        FileHeaderPKD1 h = read_header(hs);
        h.reserved1 += add;
    
        // fixed offset because layout is exactly:
        // 10 (gzip) + 2 (XLEN) + 2 (SI) + 2 (LEN)
        hs.seekp(10 + 2 + 2 + 2, std::ios::beg);
        write_header(hs, h);
        hs.flush();
    }
}

static inline void deflate_raw_to_stream(std::ostream& os,
                                         const void* src,
                                         std::size_t src_len,
                                         int level)
{
    // zlib's z_stream::avail_in is uInt (usually 32-bit). Do NOT truncate silently.
    if (src_len > static_cast<std::size_t>(std::numeric_limits<uInt>::max()))
        throw std::runtime_error("deflate_raw_to_stream: input too large for zlib ( > 4GiB )");

    z_stream zs{};
    if (deflateInit2(&zs, level, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit2 raw failed");

    zs.next_in  = reinterpret_cast<Bytef*>(const_cast<void*>(src));
    zs.avail_in = static_cast<uInt>(src_len);

    unsigned char outbuf[64 * 1024];
    int rc;
    do
    {
        zs.next_out  = outbuf;
        zs.avail_out = static_cast<uInt>(sizeof(outbuf));

        rc = deflate(&zs, zs.avail_in ? Z_NO_FLUSH : Z_FINISH);
        if (rc != Z_OK && rc != Z_STREAM_END)
        {
            deflateEnd(&zs);
            throw std::runtime_error("deflate error");
        }

        const std::size_t produced = sizeof(outbuf) - zs.avail_out;
        if (produced)
            os.write(reinterpret_cast<const char*>(outbuf),
                     static_cast<std::streamsize>(produced));
    }
    while (rc != Z_STREAM_END);

    deflateEnd(&zs);
}

inline void write_records_gzip(const std::string& path,
                               const std::vector<PackedParticle>& rec,
                               std::uint64_t released_count,
                               int level = 6,
                               bool append = false)
{
    const std::size_t bytes = rec.size() * sizeof(PackedParticle);

    // zlib crc32() length parameter is uInt -> avoid silent truncation
    if (bytes > static_cast<std::size_t>(std::numeric_limits<uInt>::max()))
        throw std::runtime_error("write_records_gzip: buffer too large for zlib ( > 4GiB )");

    if (!append)
    {
        std::ofstream os(path, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!os)
            throw std::runtime_error("write_records_gzip: cannot open " + path);

        FileHeaderPKD1 h{};
        h.reserved0 = released_count;
        h.reserved1 = static_cast<uint64_t>(rec.size());

        gzdetail::write_gzip_header(os, h, /*with_extra=*/true);

        // CRC32 over uncompressed data
        uLong crc = crc32(0L, Z_NULL, 0);
        crc = crc32(crc,
                    reinterpret_cast<const Bytef*>(rec.data()),
                    static_cast<uInt>(bytes));

        // gzip trailer ISIZE is modulo 2^32
        const uint32_t isize32 = static_cast<uint32_t>(static_cast<uint64_t>(bytes));

        // raw deflate payload
        deflate_raw_to_stream(os, rec.data(), bytes, level);
        if (!os.good())
            throw std::runtime_error("write_records_gzip: ostream error during deflate");

        // gzip trailer: CRC32 + ISIZE
        gzdetail::put_le32(os, static_cast<uint32_t>(crc));
        gzdetail::put_le32(os, isize32);
        return;
    }

    // append path
    if (rec.empty())
        return;

    // bump record count in first gzip member
    gzdetail::bump_gzip_extra_reserved1(path, static_cast<uint64_t>(rec.size()));

    std::ofstream os(path, std::ios::binary | std::ios::out | std::ios::app);
    if (!os)
        throw std::runtime_error("write_records_gzip: cannot open for append " + path);

    // appended member has no extra
    gzdetail::write_gzip_header(os, FileHeaderPKD1{}, /*with_extra=*/false);

    uLong crc = crc32(0L, Z_NULL, 0);
    crc = crc32(crc,
                reinterpret_cast<const Bytef*>(rec.data()),
                static_cast<uInt>(bytes));

    const uint32_t isize32 = static_cast<uint32_t>(static_cast<uint64_t>(bytes));

    deflate_raw_to_stream(os, rec.data(), bytes, level);
    if (!os.good())
        throw std::runtime_error("write_records_gzip: ostream error during deflate");

    gzdetail::put_le32(os, static_cast<uint32_t>(crc));
    gzdetail::put_le32(os, isize32);
}

inline std::vector<PackedParticle> read_records_gzip(const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("read_records_gzip: cannot open " + path);

    // Rewind-safe peek for PKD extra (optional)
    const std::streampos start = is.tellg();
    unsigned char hdr[10]{};
    if (!is.read(reinterpret_cast<char*>(hdr), 10))
        throw std::runtime_error("read_records_gzip: truncated header");

    if (hdr[0]==0x1f && hdr[1]==0x8b && hdr[2]==8)
    {
        const uint8_t flg = hdr[3];
        if (flg & 0x04) // FEXTRA present
        {
            const uint16_t xlen = gzdetail::get_le16(is);
            if (xlen < 4)
                throw std::runtime_error("read_records_gzip: corrupt FEXTRA (xlen < 4)");

            // We only inspect the *first* subfield. If it's not ours, we skip all extras.
            char si1 = 0, si2 = 0;
            is.get(si1); is.get(si2);
            const uint16_t len = gzdetail::get_le16(is);

            if (si1=='P' && si2=='K' && len >= sizeof(FileHeaderPKD1))
            {
                // Read our header payload
                FileHeaderPKD1 h{};
                read_all(is, &h);

                if (h.magic[0] != 'P' || h.magic[1] != 'K' || h.magic[2] != 'D')
                    throw std::runtime_error("read_records_gzip: bad PKD header");

                // Skip any extra bytes inside this subfield
                const std::streamoff extra_in_subfield =
                    static_cast<std::streamoff>(len - sizeof(FileHeaderPKD1));
                if (extra_in_subfield > 0)
                    is.seekg(extra_in_subfield, std::ios::cur);

                // Skip remaining extra subfields if present
                const uint16_t consumed = static_cast<uint16_t>(2 + 2 + len); // SI1+SI2+LEN+payload
                if (xlen > consumed)
                    is.seekg(static_cast<std::streamoff>(xlen - consumed), std::ios::cur);
            }
            else
            {
                // Not our subfield -> skip the rest of FEXTRA.
                // We already consumed 4 bytes (SI1,SI2,LEN). Need to skip remaining (xlen - 4).
                is.seekg(static_cast<std::streamoff>(xlen - 4), std::ios::cur);
            }
        }
    }

    // rewind so zlib sees the full gzip stream from the beginning
    is.clear();
    is.seekg(start);

    z_stream zs{};
    if (inflateInit2(&zs, 15 + 32) != Z_OK)
        throw std::runtime_error("inflateInit2 failed");

    std::vector<char> inbuf(64*1024), outbuf(64*1024), raw;
    zs.next_in = Z_NULL;
    zs.avail_in = 0;

    for (;;)
    {
        if (zs.avail_in == 0)
        {
            is.read(inbuf.data(), static_cast<std::streamsize>(inbuf.size()));
            const std::streamsize got = is.gcount();
            if (got == 0) break;

            zs.next_in  = reinterpret_cast<Bytef*>(inbuf.data());
            zs.avail_in = static_cast<uInt>(got);
        }

        zs.next_out  = reinterpret_cast<Bytef*>(outbuf.data());
        zs.avail_out = static_cast<uInt>(outbuf.size());

        const int ret = inflate(&zs, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END)
        {
            inflateEnd(&zs);
            throw std::runtime_error("inflate error: " + std::to_string(ret));
        }

        const std::size_t produced = outbuf.size() - zs.avail_out;
        if (produced)
            raw.insert(raw.end(), outbuf.data(), outbuf.data() + produced);

        if (ret == Z_STREAM_END)
        {
            // More concatenated members?
            if (zs.avail_in > 0 || is.peek() != EOF)
            {
                if (inflateReset2(&zs, 15 + 32) != Z_OK)
                {
                    inflateEnd(&zs);
                    throw std::runtime_error("inflateReset2 failed");
                }
                continue;
            }
            break;
        }
    }

    inflateEnd(&zs);

    if (raw.size() % sizeof(PackedParticle) != 0)
        throw std::runtime_error("read_records_gzip: size mismatch");

    std::vector<PackedParticle> rec(raw.size() / sizeof(PackedParticle));
    std::memcpy(rec.data(), raw.data(), raw.size());
    return rec;
}

#endif // PK_IO_USE_ZLIB
