// fixed_point_core.hpp
// Core math + quantization + packed structs + encode/decode.
#pragma once
#include <cstdint>
#include <cmath>
#include <array>
#include <limits>
#include <numbers>      // std::numbers::pi_v

// ------------------------ constants ------------------------
// Keep a single source of truth for scaling across PKD1/PKD2 and seed loaders.
inline constexpr double   fp_pi             = std::numbers::pi_v<double>;
inline constexpr uint32_t fp_ticks_per_rad  = 1'000'000u;     // ticks per radian
inline constexpr double   dbl_ticks_per_rad  = 1'000'000.;    // ticks per radian (double)
inline constexpr uint64_t fp_ticks_per_sec   = 1u;             // seconds per tick (PKD1)

// ------------------------ angle utils ------------------------
inline double wrap_lon(double r)
{
    r = std::fmod(r,2.*fp_pi);
    if (r<=-fp_pi) r+=2.*fp_pi;
    else if (r>fp_pi) r-=2.*fp_pi;
    return r;
}

// ------------------------ quantization -----------------------
// Position: ticks = radians * fp_ticks_per_rad
template <typename Real>
inline int32_t quant_coord(Real rad)
{
    long double v = std::llround(static_cast<long double>(rad) * static_cast<long double>(fp_ticks_per_rad));
    if (v>std::numeric_limits<int32_t>::max()) v = std::numeric_limits<int32_t>::max();
    if (v<std::numeric_limits<int32_t>::min()) v = std::numeric_limits<int32_t>::min();
    return static_cast<int32_t>(v);
}

template <typename Real>
inline Real dequant_coord(int32_t q)
{
    return static_cast<Real>(q) / static_cast<Real>(fp_ticks_per_rad);
}

// Time: 1 tick = 1 second by default (PKD1 helpers)
template <typename Real>
inline int64_t quant_time(Real seconds)
{
    return static_cast<int64_t>(std::llround(static_cast<long double>(seconds)));
}

template <typename Real>
inline Real dequant_time(int64_t ticks)
{
    return static_cast<Real>(ticks);
}

// ---------------------- packed record (PKD1) ------------------------
#pragma pack(push,1)
struct PackedParticle
{
    int32_t lon0_q;  // ticks (micro-rad)
    int32_t lat0_q;  // ticks
    int32_t lon1_q;  // ticks
    int32_t lat1_q;  // ticks
    int64_t t0_s;    // seconds since reference
    int64_t t1_s;    // seconds since reference
    uint8_t flags;
};
#pragma pack(pop)

// ------------------- PKD1 file header (struct only) ------------------
struct FileHeaderPKD1
{
    std::array<char,8> magic = {'P','K','D','1',0,0,0,0};
    uint32_t version   = 1;
    uint32_t pos_scale = fp_ticks_per_rad;  // ticks per rad
    uint64_t time_tick = fp_ticks_per_sec;  // seconds per tick
    uint64_t reserved0 = 0;
    uint64_t reserved1 = 0;
};

// ---------------- encode/decode one record -----------------
template<typename PosT,typename TimeT>
inline PackedParticle encode_record(PosT lon0, PosT lat0, TimeT t0,
                                    PosT lon1, PosT lat1, TimeT t1,
                                    uint8_t flags = 0)
{
    PackedParticle r{};
    r.lon0_q = quant_coord(wrap_lon(lon0));
    r.lat0_q = quant_coord(lat0);
    r.lon1_q = quant_coord(wrap_lon(lon1));
    r.lat1_q = quant_coord(lat1);
    r.t0_s   = quant_time(t0);
    r.t1_s   = quant_time(t1);
    r.flags  = flags;
    return r;
}

template<typename PosT,typename TimeT>
inline void decode_record(const PackedParticle& r,
                          PosT& lon0, PosT& lat0, TimeT& t0,
                          PosT& lon1, PosT& lat1, TimeT& t1,
                          uint8_t& flags)
{
    lon0 = dequant_coord<PosT>(r.lon0_q);
    lat0 = dequant_coord<PosT>(r.lat0_q);
    lon1 = dequant_coord<PosT>(r.lon1_q);
    lat1 = dequant_coord<PosT>(r.lat1_q);
    t0   = dequant_time<TimeT>(r.t0_s);
    t1   = dequant_time<TimeT>(r.t1_s);
    flags = r.flags;
}
