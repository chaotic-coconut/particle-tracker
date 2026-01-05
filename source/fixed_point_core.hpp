/********************************************************************
 *  fixed_point_core.hpp  –  basic packing helpers (no compression) *
 ********************************************************************/
#pragma once
#include <cstdint>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstring>   // std::memcpy
#include <algorithm> // std::min

// ------------------------ angle utils ------------------------
inline constexpr double pi=3.141592653589793238462643383279502884;

inline double wrap_lon(double r)
{
	r=std::fmod(r,2*pi);
	if (r<=-pi) r+=2*pi;
	else if (r>pi) r-=2*pi;
	return r;
}

// ------------------------ quantization -----------------------
// Position scale: micro-radians (approx. 6.37 m per least significant bit on Earth)
template <typename DataType>
inline int32_t quant_coord(DataType rad)
{
	return static_cast<int32_t>(std::llround(static_cast<long double>(rad)*1'000'000.L));
}

template <typename DataType>
inline DataType dequant_coord(int32_t q)
{
	return static_cast<DataType>(q)/static_cast<DataType>(1'000'000.);
}

// Time scale: 1 tick = 1 second
template <typename TimeType>
inline int64_t quant_time(TimeType seconds)
{
	return static_cast<int64_t>(std::llround(static_cast<long double>(seconds)));
}

template <typename TimeType>
inline TimeType dequant_time(int64_t ticks)
{
	return static_cast<TimeType>(ticks);
}

// ---------------------- packed record ------------------------
#pragma pack(push,1)
struct PackedParticle
{
	int32_t lon0_q;  // micro-rad
	int32_t lat0_q;  // micro-rad
	int32_t lon1_q;  // micro-rad
	int32_t lat1_q;  // micro-rad
	int64_t t0_s;    // seconds since reference
	int64_t t1_s;    // seconds since reference
	uint8_t flags;   //
};
#pragma pack(pop)

// -------------------- simple file header ---------------------
struct FileHeader
{
	std::array<char,8> magic={'P','K','D','1',0,0,0,0 }; // “PKD1\0\0\0\0”
	uint32_t version  =1;
	uint32_t pos_scale=1'000'000;  // µrad per LSB
	uint64_t time_tick=1;          // s  per tick
	uint64_t reserved0=0;
	uint64_t reserved1=0;
};

// ------------------- safe read/write bytes -------------------
template <class T>
inline void write_all(std::ostream& os,const T* ptr,std::size_t count=1)
{
	os.write(reinterpret_cast<const char*>(ptr),sizeof(T)*count);
	if (!os) throw std::runtime_error("write_all: write failed");
}

template <class T>
inline void read_all(std::istream& is,T* ptr,std::size_t count=1)
{
	is.read(reinterpret_cast<char*>(ptr),sizeof(T)*count);
	if (!is) throw std::runtime_error("read_all: read failed");
}


inline void write_header(std::ostream& os,const FileHeader& h)
{
	write_all(os,&h);
}

inline FileHeader read_header(std::istream& is)
{
	FileHeader h{};
	read_all(is,&h);
	if (h.magic[0]!='P' || h.magic[1]!='K' || h.magic[2]!='D')
		throw std::runtime_error("read_header: bad magic");
	if (h.version!=1)
		throw std::runtime_error("read_header: unsupported version");
	return h;
}

// ---------------- encode/decode one record -----------------
template<typename PosT,typename TimeT>
inline PackedParticle encode_record(PosT lon0,PosT lat0,TimeT t0,
                                    PosT lon1,PosT lat1,TimeT t1,
                                    uint8_t flags=0)
{
	PackedParticle r{};
	r.lon0_q=quant_coord(wrap_lon(lon0));
	r.lat0_q=quant_coord(lat0);
	r.lon1_q=quant_coord(wrap_lon(lon1));
	r.lat1_q=quant_coord(lat1);
	r.t0_s  =quant_time(t0);
	r.t1_s  =quant_time(t1);
	r.flags =flags;
	return r;
}

template<typename PosT,typename TimeT>
inline void decode_record(const PackedParticle& r,
                          PosT& lon0,PosT& lat0,TimeT& t0,
                          PosT& lon1,PosT& lat1,TimeT& t1,
                          uint8_t& flags)
{
	lon0=dequant_coord<PosT>(r.lon0_q);
	lat0=dequant_coord<PosT>(r.lat0_q);
	lon1=dequant_coord<PosT>(r.lon1_q);
	lat1=dequant_coord<PosT>(r.lat1_q);
	t0  =dequant_time<TimeT>(r.t0_s);
	t1  =dequant_time<TimeT>(r.t1_s);
	flags=r.flags;
}

/* ---------------- plain (uncompressed) I/O ----------------- */
inline void write_records_binary(const std::string& path,
                                 const std::vector<PackedParticle>& rec,
                                 bool append=false)
{
	std::ios::openmode mode=std::ios::binary | std::ios::out;
	if (append) mode |=std::ios::app;

	std::ofstream os(path,mode);
	if (!os) throw std::runtime_error("cannot open "+path);

	if (!append) {FileHeader h;write_header(os,h);}

	if (!rec.empty()) write_all(os,rec.data(),rec.size());
}

inline std::vector<PackedParticle>
read_records_binary(const std::string& path)
{
	std::ifstream is(path,std::ios::binary);
	if (!is) throw std::runtime_error("read_records_binary: cannot open file "+path);

	(void)read_header(is);

	is.seekg(0,std::ios::end);
	std::streamoff end=is.tellg();
	is.seekg(sizeof(FileHeader),std::ios::beg);
	std::streamoff bytes=end-std::streamoff(sizeof(FileHeader));

	if (bytes<0 || bytes%static_cast<std::streamoff>(sizeof(PackedParticle))!=0)
		throw std::runtime_error("file size mismatch");

	std::size_t n=static_cast<std::size_t>(bytes/sizeof(PackedParticle));
	std::vector<PackedParticle> rec(n);
	if (n) read_all(is,rec.data(),n);
	return rec;
}

