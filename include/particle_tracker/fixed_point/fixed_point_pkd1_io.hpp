// fixed_point_pkd1_io.hpp
// PKD1 plain (uncompressed) file I/O kept out of fixed_point_core.hpp.
#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>
#include "fixed_point_core.hpp" // FileHeaderPKD1, PackedParticle

// ------------------- safe read/write bytes -------------------
template <class T>
inline void write_all(std::ostream& os, const T* ptr, std::size_t count = 1)
{
    os.write(reinterpret_cast<const char*>(ptr), sizeof(T)*count);
    if (!os) throw std::runtime_error("write_all: write failed");
}

template <class T>
inline void read_all(std::istream& is, T* ptr, std::size_t count = 1)
{
    is.read(reinterpret_cast<char*>(ptr), sizeof(T) * count);
    if (!is) throw std::runtime_error("read_all: read failed");
}

inline void write_header(std::ostream& os, const FileHeaderPKD1& h)
{
    write_all(os, &h);
}

inline FileHeaderPKD1 read_header(std::istream& is)
{
    FileHeaderPKD1 h{};
    read_all(is, &h);

    if (h.magic[0] != 'P' || h.magic[1] != 'K' || h.magic[2] != 'D' || h.magic[3] != '1')
        throw std::runtime_error("read_header: bad magic (expected PKD1)");
    if (h.version != 1)
        throw std::runtime_error("read_header: unsupported PKD1 version");

    return h;
}

/* ---------------- plain (uncompressed) I/O ----------------- */
inline void write_records_binary(const std::string& path,
                                 const std::vector<PackedParticle>& rec,
                                 bool append = false)
{
    std::ios::openmode mode = std::ios::binary | std::ios::out;
    if (append) mode |= std::ios::app;

    std::ofstream os(path, mode);
    if (!os) throw std::runtime_error("cannot open " + path);

    if (!append)
    {
        FileHeaderPKD1 h{};
        write_header(os, h);
    }

    if (!rec.empty()) write_all(os, rec.data(), rec.size());
}

inline std::vector<PackedParticle> read_records_binary(const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("read_records_binary: cannot open file " + path);

    (void)read_header(is);

    is.seekg(0, std::ios::end);
    std::streamoff end = is.tellg();
    is.seekg(static_cast<std::streamoff>(sizeof(FileHeaderPKD1)), std::ios::beg);
    std::streamoff bytes = end - static_cast<std::streamoff>(sizeof(FileHeaderPKD1));

    if (bytes < 0 || bytes % static_cast<std::streamoff>(sizeof(PackedParticle)) != 0)
        throw std::runtime_error("read_records_binary: file size mismatch");

    std::size_t n = static_cast<std::size_t>(bytes / sizeof(PackedParticle));
    std::vector<PackedParticle> rec(n);
    if (n) read_all(is, rec.data(), n);
    return rec;
}

