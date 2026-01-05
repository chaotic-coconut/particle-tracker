/********************************************************************
 *  fixed_point_io_codec.hpp  –  compressed I/O helpers             *
 *                                                                   *
 *  Requires   fixed_point_core.hpp                                  *
 *  Define ONE of:                                                   *
 *        PK_IO_USE_ZSTD   – use Zstandard                           *
 *        PK_IO_USE_ZLIB   – use zlib / gzip                         *
 ********************************************************************/
#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <cstring>	// std::memcpy

/*------------------------------------------------------------------*/
/*  If the build system forgets the macro we abort at compile time  */
/*------------------------------------------------------------------*/
#if !defined(PK_IO_USE_ZSTD) && !defined(PK_IO_USE_ZLIB)
#  error "Define either PK_IO_USE_ZSTD or PK_IO_USE_ZLIB before including fixed_point_io_codec.hpp"
#endif
#if defined(PK_IO_USE_ZSTD) && defined(PK_IO_USE_ZLIB)
#  error "Define only ONE of PK_IO_USE_ZSTD or PK_IO_USE_ZLIB, not both"
#endif

/*------------------------------------------------------------------*/
/*  Zstandard implementation                                        */
/*------------------------------------------------------------------*/
#ifdef PK_IO_USE_ZSTD
#include <zstd.h>     // link with -lzstd

// --- zstd writer that stores FileHeader in a skippable frame ---
inline void write_records_zstd(const std::string& path,
                               const std::vector<PackedParticle>& rec,
                               std::uint64_t released_count,
                               int level=3,
                               bool append=false)
{
	// 1) On first write, create skippable frame carrying FileHeader
	if (!append)
	{
		std::ofstream os(path,std::ios::binary | std::ios::out | std::ios::trunc);
		if (!os) throw std::runtime_error("write_records_zstd: cannot open "+path);
		FileHeader h{};
		h.reserved0=released_count;                  // total released that day
		h.reserved1=static_cast<uint64_t>(rec.size()); // first batch stored

		// skippable frame magic (choose any 0x184D2A50..57)
		const uint32_t magic=0x184D2A50u;
		const uint32_t size =static_cast<uint32_t>(sizeof(FileHeader));
		os.write(reinterpret_cast<const char*>(&magic),4);
		os.write(reinterpret_cast<const char*>(&size), 4);
		write_header(os, h);
		os.flush();
	}
	else
	{
		// bump reserved1 inside skippable frame at file start
		std::fstream hs(path,std::ios::binary | std::ios::in | std::ios::out);
		if (!hs) throw std::runtime_error("write_records_zstd: cannot reopen for header "+path);
		uint32_t magic=0,size=0;
		hs.read(reinterpret_cast<char*>(&magic),4);
		hs.read(reinterpret_cast<char*>(&size), 4);
		if ((magic & 0xFFFFFFF0u)!=0x184D2A50u || size<sizeof(FileHeader))
			throw std::runtime_error("write_records_zstd: no skippable header found");
		FileHeader h=read_header(hs);
		h.reserved1+=static_cast<uint64_t>(rec.size());
		hs.seekp(8,std::ios::beg);          // after magic+size
		write_header(hs,h);
		hs.flush();
	}

	if (rec.empty()) return;

	// 2) Append a normal zstd frame with the records
	std::ofstream os(path,std::ios::binary | std::ios::out | std::ios::app);
	if (!os) throw std::runtime_error("write_records_zstd: cannot open for append "+path);

	const size_t src_bytes=rec.size()*sizeof(PackedParticle);
	const size_t bound    =ZSTD_compressBound(src_bytes);
	std::vector<char> cbuf(bound);

	const size_t rc=ZSTD_compress(cbuf.data(),bound,rec.data(),src_bytes,level);
	if (ZSTD_isError(rc))
		throw std::runtime_error("write_records_zstd: "+std::string(ZSTD_getErrorName(rc)));

	write_all(os,cbuf.data(), rc);
}

inline std::vector<PackedParticle>
read_records_zstd(const std::string& path)
{
	std::ifstream is(path,std::ios::binary);
	if (!is) throw std::runtime_error("read_records_zstd: cannot open "+path);

	// Try to read a skippable frame carrying FileHeader
	std::streampos start=is.tellg();
	uint32_t magic=0,size=0;
	is.read(reinterpret_cast<char*>(&magic),4);
	is.read(reinterpret_cast<char*>(&size), 4);

	if (is && ((magic & 0xFFFFFFF0u) == 0x184D2A50u))
	{
		if (size<sizeof(FileHeader))
			throw std::runtime_error("read_records_zstd: skippable too small");
		FileHeader h{};
		read_all(is,&h);
		if (h.magic[0]!='P' || h.magic[1]!='K' || h.magic[2]!='D')
			throw std::runtime_error("read_records_zstd: bad PKD header");
		if (size>sizeof(FileHeader))
			is.seekg(static_cast<std::streamoff>(size-sizeof(FileHeader)),std::ios::cur);
		// Now positioned at the first zstd frame
	}
	else
	{
		// Backward compatibility: old files had a raw FileHeader at the start.
		is.clear();
		is.seekg(start);
		(void)read_header(is);
	}

	ZSTD_DCtx* dctx=ZSTD_createDCtx();
	if (!dctx) throw std::runtime_error("read_records_zstd: createDCtx failed");

	const size_t in_chunk =ZSTD_DStreamInSize();
	const size_t out_chunk=ZSTD_DStreamOutSize();
	std::vector<char> inbuf(in_chunk),outbuf(out_chunk),raw;

	ZSTD_inBuffer  in  {inbuf.data(),0,0};
	ZSTD_outBuffer out {outbuf.data(),outbuf.size(),0};

	while (true)
	{
		if (in.pos==in.size)
		{
		is.read(inbuf.data(),static_cast<std::streamsize>(inbuf.size()));
		in.size=static_cast<size_t>(is.gcount());
		in.pos =0;
		if (in.size==0) break; // EOF
		}
		out.pos=0;
		size_t ret=ZSTD_decompressStream(dctx,&out,&in);
		if (ZSTD_isError(ret))
		{
			ZSTD_freeDCtx(dctx);
			throw std::runtime_error("read_records_zstd: "+std::string(ZSTD_getErrorName(ret)));
		}
		if (out.pos) raw.insert(raw.end(),outbuf.data(),outbuf.data()+out.pos);
	}
	ZSTD_freeDCtx(dctx);

	if (raw.size() % sizeof(PackedParticle)!= 0)
		throw std::runtime_error("read_records_zstd: size mismatch");

	std::vector<PackedParticle> rec(raw.size()/sizeof(PackedParticle));
	std::memcpy(rec.data(),raw.data(),raw.size());
	return rec;
}

#endif  // PK_IO_USE_ZSTD


/*------------------------------------------------------------------*/
/*  zlib / gzip implementation                                      */
/*------------------------------------------------------------------*/
#ifdef PK_IO_USE_ZLIB
#include <zlib.h>	// link with -lz

/* helper – throw if zlib error */
inline void zerr(int code,const char* where)
{
	if (code!=Z_OK && code!=Z_STREAM_END)
		throw std::runtime_error(std::string(where)+" : zlib error "+std::to_string(code));
}

namespace gzdetail
{
	inline void put_le16(std::ostream& os,uint16_t v)
	{
		unsigned char b[2]={(unsigned char)(v & 0xFFu),(unsigned char)((v>>8) & 0xFFu) };
		os.write(reinterpret_cast<const char*>(b),2);
	}
	inline void put_le32(std::ostream& os,uint32_t v)
	{
		unsigned char b[4]={(unsigned char)(v & 0xFFu),
				    (unsigned char)((v>>8) & 0xFFu),
				    (unsigned char)((v>>16)& 0xFFu),
				    (unsigned char)((v>>24)& 0xFFu)};
		os.write(reinterpret_cast<const char*>(b),4);
	}

	// Write a gzip header with FEXTRA that contains a PKD subfield with FileHeader bytes
	inline void write_gzip_header(std::ostream& os,const FileHeader& h,bool with_extra=true)
	{
		const uint8_t ID1=0x1f,ID2=0x8b,CM=8;
		uint8_t FLG=with_extra ? 0x04 : 0x00;  // FEXTRA
		os.put(ID1);os.put(ID2);os.put(CM);os.put(FLG);
		put_le32(os,0);     // MTIME
		os.put(0);           // XFL
		//os.put(255);         // OS = unknown
		unsigned char os_byte=0xFF;
		os.write(reinterpret_cast<const char*>(&os_byte),1);

		if (with_extra)
		{
			// Extra: one subfield with SI1='P', SI2='K'
			const uint16_t sublen=static_cast<uint16_t>(sizeof(FileHeader));
			const uint16_t xlen  =static_cast<uint16_t>(2+2+sublen);
			put_le16(os,xlen);
			os.put('P');os.put('K');            // SI1, SI2 (choose your own)
			put_le16(os,sublen);                // LEN
			write_header(os,h);                 // payload
		}
	}

	// Update reserved1 in the extra subfield at the start of the file
	inline void bump_gzip_extra_reserved1(const std::string& path,uint64_t add)
	{
		std::fstream hs(path,std::ios::binary | std::ios::in | std::ios::out);
		if (!hs) throw std::runtime_error("write_records_gzip: cannot reopen "+path);

		unsigned char hdr[10];
		hs.read(reinterpret_cast<char*>(hdr),10);
		if (hs.gcount()!=10 || hdr[0]!=0x1f || hdr[1]!=0x8b || hdr[2]!=8)
		throw std::runtime_error("write_records_gzip: not a gzip file");
		uint8_t flg=hdr[3];
		if ((flg & 0x04)==0) return; // no FEXTRA -> nothing to bump

		uint16_t xlen;hs.read(reinterpret_cast<char*>(&xlen), 2);

		// Read first subfield header (we wrote exactly one)
		char si1,si2;uint16_t len;
		hs.get(si1);hs.get(si2);
		hs.read(reinterpret_cast<char*>(&len),2);
		if (si1!='P' || si2!='K' || len<sizeof(FileHeader))
			throw std::runtime_error("write_records_gzip: PKD extra not found");

		// Read-modify-write FileHeader payload
		FileHeader h=read_header(hs);
		h.reserved1+=add;
		hs.seekp(10+2+2+2,std::ios::beg); // after ID..OS, XLEN, SI, LEN
		write_header(hs,h);
		hs.flush();
	}
}

// ---- helper: raw DEFLATE a buffer to a stream (no zlib/gzip wrapper) ----
static inline void deflate_raw_to_stream(std::ostream& os,const void* src,std::size_t src_len,int level)
{
	z_stream zs{};
	if (deflateInit2(&zs,level,Z_DEFLATED,-15,9,Z_DEFAULT_STRATEGY)!=Z_OK)
		throw std::runtime_error("deflateInit2 raw failed");

	zs.next_in =reinterpret_cast<Bytef*>(const_cast<void*>(src));
	zs.avail_in=static_cast<uInt>(src_len);

	unsigned char outbuf[64*1024];
	int rc;
	do
	{
		zs.next_out =outbuf;
		zs.avail_out=sizeof(outbuf);

		rc=deflate(&zs,zs.avail_in ? Z_NO_FLUSH : Z_FINISH);
		if (rc!=Z_OK && rc!=Z_STREAM_END)
		{
			deflateEnd(&zs);
			throw std::runtime_error("deflate error");
		}

		const std::size_t produced=sizeof(outbuf)-zs.avail_out;
		if (produced)
		{
			os.write(reinterpret_cast<const char*>(outbuf),
			static_cast<std::streamsize>(produced));
		}
	}
	while (rc!=Z_STREAM_END);

	deflateEnd(&zs);
}

// --- gzip writer that keeps header in FEXTRA and appends members ---
inline void write_records_gzip(const std::string& path,
                               const std::vector<PackedParticle>& rec,
                               std::uint64_t released_count,
                               int level=6,
                               bool append=false)
{
	std::cerr<<"write_records_gzip: rec="<<rec.size()<<" append="<<append<<"\n";

	if (!append)
	{
		// First member: write header with FEXTRA carrying FileHeader
		std::ofstream os(path,std::ios::binary | std::ios::out | std::ios::trunc);
		if (!os) throw std::runtime_error("write_records_gzip: cannot open "+path);
		FileHeader h{};h.reserved0=released_count;h.reserved1=static_cast<uint64_t>(rec.size());
		gzdetail::write_gzip_header(os,h,/*with_extra=*/true);

		// Deflate raw stream + trailer
		uLong crc=crc32(0L,Z_NULL,0);
		crc=crc32(crc,reinterpret_cast<const Bytef*>(rec.data()),rec.size()*sizeof(PackedParticle));
		uLong isize=static_cast<uLong>(rec.size()*sizeof(PackedParticle));

		deflate_raw_to_stream(os,rec.data(),rec.size()*sizeof(PackedParticle),level);
		if (!os.good()) throw std::runtime_error("ostream error during deflate");

		// trailer: CRC32 and ISIZE (little-endian)
		gzdetail::put_le32(os,static_cast<uint32_t>(crc));
		gzdetail::put_le32(os,static_cast<uint32_t>(isize));
	}
	else
	{
		if (rec.empty()) return; // nothing to append; don't bump

		// bump count in the first member’s PKD extra
		gzdetail::bump_gzip_extra_reserved1(path,static_cast<uint64_t>(rec.size()));

		std::ofstream os(path, std::ios::binary | std::ios::out | std::ios::app);
		if (!os) throw std::runtime_error("write_records_gzip: cannot open for append "+path);

		gzdetail::write_gzip_header(os,FileHeader{}, /*with_extra=*/false);

		uLong crc=crc32(0L,Z_NULL,0);
		crc=crc32(crc,reinterpret_cast<const Bytef*>(rec.data()),rec.size()*sizeof(PackedParticle));
		uLong isize=static_cast<uLong>(rec.size()*sizeof(PackedParticle));

		deflate_raw_to_stream(os,rec.data(),rec.size()*sizeof(PackedParticle),level);
		if (!os.good()) throw std::runtime_error("ostream error during deflate");

		gzdetail::put_le32(os,static_cast<uint32_t>(crc));
		gzdetail::put_le32(os,static_cast<uint32_t>(isize));
	}
}

/* ---------- read (gzip stream) ------------- */
inline std::vector<PackedParticle>
read_records_gzip(const std::string& path)
{
	std::ifstream is(path,std::ios::binary);
	if (!is) throw std::runtime_error("read_records_gzip: cannot open "+path);

	// Peek first member header to read PKD extra, then rewind so zlib can parse it again.
	std::streampos start=is.tellg();
	unsigned char hdr[10];
	if (!is.read(reinterpret_cast<char*>(hdr),10))
		throw std::runtime_error("read_records_gzip: truncated header");

	if (hdr[0]==0x1f && hdr[1]==0x8b && hdr[2]==8)
	{
		uint8_t flg=hdr[3];
		if (flg & 0x04)
		{
			// FEXTRA
			uint16_t xlen;is.read(reinterpret_cast<char*>(&xlen),2);

			char si1,si2;uint16_t len;
			is.get(si1);is.get(si2);
			is.read(reinterpret_cast<char*>(&len),2);

			if (si1=='P' && si2=='K' && len >= sizeof(FileHeader))
			{
				FileHeader h{};
				read_all(is,&h);
				if (h.magic[0]!='P' || h.magic[1]!='K' || h.magic[2]!='D')
				throw std::runtime_error("read_records_gzip: bad PKD header");

				// Skip any leftover subfield bytes
				if (len>sizeof(FileHeader))
				is.seekg(static_cast<std::streamoff>(len-sizeof(FileHeader)),std::ios::cur);

				// Skip any other extra subfields if present
				int consumed=2+2+len; // SI+LEN+payload
				if (xlen>consumed)
					is.seekg(static_cast<std::streamoff>(xlen-consumed),std::ios::cur);
			}
			else
			{
				// No PKD subfield; skip the whole extra block
				is.seekg(static_cast<std::streamoff>(xlen-4),std::ios::cur);
			}
		}
	}
	// Rewind so zlib sees the full gzip member(s)
	is.clear();
	is.seekg(start);

	z_stream zs{};
	if (inflateInit2(&zs,15+32)!=Z_OK) // auto zlib/gzip with header & trailer
		throw std::runtime_error("inflateInit2 failed");

	std::vector<char> inbuf(64*1024),outbuf(64*1024),raw;
	zs.next_in=Z_NULL;
	zs.avail_in=0;

	for (;;)
	{
		if (zs.avail_in==0)
		{
			is.read(inbuf.data(),static_cast<std::streamsize>(inbuf.size()));
			std::streamsize got=is.gcount();
			if (got==0) break; // EOF
			zs.next_in =reinterpret_cast<Bytef*>(inbuf.data());
			zs.avail_in=static_cast<uInt>(got);
		}

		zs.next_out=reinterpret_cast<Bytef*>(outbuf.data());
		zs.avail_out=static_cast<uInt>(outbuf.size());

		int ret=inflate(&zs,Z_NO_FLUSH);
		if (ret!=Z_OK && ret!=Z_STREAM_END)
			throw std::runtime_error("inflate error: "+std::to_string(ret));

		std::size_t produced=outbuf.size()-zs.avail_out;
		if (produced)
		raw.insert(raw.end(),outbuf.data(),outbuf.data()+produced);

		if (ret==Z_STREAM_END)
		{
			// More concatenated members?
			if (zs.avail_in > 0 || is.peek() != EOF)
			{
				if (inflateReset2(&zs,15+32)!=Z_OK)
				throw std::runtime_error("inflateReset2 failed");
				continue;
			}
		break;
		}
	}
	inflateEnd(&zs);

	if (raw.size() % sizeof(PackedParticle)!=0)
		throw std::runtime_error("read_records_gzip: size mismatch");

	std::vector<PackedParticle> rec(raw.size()/sizeof(PackedParticle));
	std::memcpy(rec.data(),raw.data(),raw.size());
	return rec;
}
#endif  // PK_IO_USE_ZLIB

