#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <cstdio>
#include <limits>

#include <zlib.h>

#ifdef PK_IO_USE_ZSTD
#  error "PK_IO_USE_ZSTD is unsupported; PKD2 files use zlib compression."
#endif
#ifndef PK_IO_USE_ZLIB
#  error "Define PK_IO_USE_ZLIB before including seed_loader_v2.hpp."
#endif
#include "particle_tracker/fixed_point/fixed_point_core.hpp"   // wrap_lon(), pi
#include "particle_tracker/fixed_point/fixed_point_pkd2.hpp"   // pkd2::{Header,TocEntry,BlockHeader}

// =======================
// Public API
// =======================

// MUST match Writer(pos_scale=1'000'000u) used for PKD2 trajectories.
inline constexpr double POS_TICKS_PER_RAD = 1'000'000.0;

enum class SeedPoint { First, Last };

// Seed record:
// - ACTIVE seed point: lon_rad/lat_rad/t_start_s  (this is what simulation starts from)
// - STOP target (optional): lon_stop_rad/lat_stop_rad/t_stop_s
// - For convenience/back-compat: lon_start_rad/lat_start_rad mirror ACTIVE.
struct SeedInit
{
    uint64_t seed_id{};

    // ACTIVE point
    double   lon_rad{};    // radians in (-pi,pi]
    double   lat_rad{};    // radians
    int64_t  t_start_s{};  // seconds (t2000 for PKD2; CSV depends on caller)

    // Convenience/back-compat: mirrors ACTIVE
    double   lon_start_rad{};
    double   lat_start_rad{};

    // Optional STOP target
    double   lon_stop_rad{};
    double   lat_stop_rad{};
    int64_t  t_stop_s{};
    bool     has_stop{false};
};

using SeedRec = SeedInit;

// =======================
// Internal helpers
// =======================
namespace seed_v2_detail
{
    inline bool ends_with_ci(const std::string& s, const std::string& suf)
    {
        if (s.size() < suf.size()) return false;
        auto lower = [](unsigned char c){ return static_cast<char>(std::tolower(c)); };

        std::string a = s.substr(s.size() - suf.size());
        std::transform(a.begin(), a.end(), a.begin(), lower);

        std::string b = suf;
        std::transform(b.begin(), b.end(), b.begin(), lower);

        return a == b;
    }

    inline void fread_or_throw(FILE* f, void* p, size_t n, const char* what)
    {
        if (std::fread(p, 1, n, f) != n)
            throw std::runtime_error(std::string("PKD2: fread failed: ") + what);
    }

    inline void fseek_or_throw(FILE* f, long off, int whence, const char* what)
    {
        if (std::fseek(f, off, whence) != 0)
            throw std::runtime_error(std::string("PKD2: fseek failed: ") + what);
    }

    inline std::vector<unsigned char> zlib_decompress_exact(const unsigned char* comp,
                                                            size_t comp_sz,
                                                            size_t raw_sz)
    {
        std::vector<unsigned char> raw(raw_sz);
        uLongf out = static_cast<uLongf>(raw_sz);
        int rc = ::uncompress(raw.data(), &out, comp, static_cast<uLong>(comp_sz));
        if (rc != Z_OK || out != raw_sz)
            throw std::runtime_error("PKD2: zlib uncompress failed (wrong codec or corrupt block)");
        return raw;
    }

    // Decode first+last sample of one PKD2 block.
    // Returns two "ACTIVE-only" records: first and last (no stop metadata assigned here).
    inline void decode_block_first_last(const pkd2::Header& hdr,
                                        const pkd2::TocEntry& te,
                                        FILE* f,
                                        SeedInit& first_out,
                                        SeedInit& last_out)
    {
        std::vector<unsigned char> comp(te.comp_size);
        fseek_or_throw(f, static_cast<long>(te.file_offset), SEEK_SET, "block");
        fread_or_throw(f, comp.data(), comp.size(), "block bytes");

        auto raw = zlib_decompress_exact(comp.data(), comp.size(), te.raw_size);
        if (raw.size() < sizeof(pkd2::BlockHeader))
            throw std::runtime_error("PKD2: raw block too small");

        const auto* bh = reinterpret_cast<const pkd2::BlockHeader*>(raw.data());
        const uint32_t n = bh->n_points;
        const uint16_t stride = bh->stride ? bh->stride : 1;
        if (n == 0) throw std::runtime_error("PKD2: n_points==0 (invalid)");

        // Strong sanity: pos_scale should match your writer
        if (hdr.pos_scale_ticks_per_rad != static_cast<uint32_t>(POS_TICKS_PER_RAD))
            throw std::runtime_error("PKD2: pos_scale_ticks_per_rad mismatch (file vs expected 1e6)");

        const double inv = 1.0 / static_cast<double>(hdr.pos_scale_ticks_per_rad);

        const double L0 = static_cast<double>(bh->lon0_q) * inv;
        const double B0 = static_cast<double>(bh->lat0_q) * inv;

        const int64_t t0 = hdr.t0_epoch_s
                         + static_cast<int64_t>(bh->start_index) * static_cast<int64_t>(hdr.dt_s);

        const int32_t* dlon = reinterpret_cast<const int32_t*>(raw.data() + sizeof(pkd2::BlockHeader));
        const int32_t* dlat = dlon + (n > 1 ? (n - 1) : 0);

        double L = L0;
        double B = B0;
        for (uint32_t i = 1; i < n; ++i)
        {
            L += static_cast<double>(dlon[i - 1]) * inv;
            B += static_cast<double>(dlat[i - 1]) * inv;
        }

        const int64_t t_last = t0
                             + static_cast<int64_t>(n - 1)
                               * static_cast<int64_t>(hdr.dt_s)
                               * static_cast<int64_t>(stride);

        // FIRST active
        first_out = {};
        first_out.seed_id   = bh->seed_id;
        first_out.lon_rad   = wrap_lon(L0);
        first_out.lat_rad   = B0;
        first_out.t_start_s = t0;
        first_out.lon_start_rad = first_out.lon_rad;
        first_out.lat_start_rad = first_out.lat_rad;
        first_out.has_stop  = false;

        // LAST active
        last_out = {};
        last_out.seed_id    = bh->seed_id;
        last_out.lon_rad    = wrap_lon(L);
        last_out.lat_rad    = B;
        last_out.t_start_s  = t_last;
        last_out.lon_start_rad = last_out.lon_rad;
        last_out.lat_start_rad = last_out.lat_rad;
        last_out.has_stop   = false;
    }

    // PKD2 loader: attach stop target depending on which_point.
    inline std::vector<SeedInit> load_seeds_from_pkd2_zlib(const std::string& path, SeedPoint which_point)
    {
        FILE* f = std::fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("PKD2: cannot open " + path);

        pkd2::Header hdr{};
        try
        {
            fread_or_throw(f, &hdr, sizeof(hdr), "header");

            if (!(hdr.magic[0]=='P' && hdr.magic[1]=='K' && hdr.magic[2]=='D' && hdr.magic[3]=='2'))
                throw std::runtime_error("PKD2: bad magic");
            if (hdr.version != 3)
                throw std::runtime_error("PKD2: unsupported version (expected 3)");
            if (hdr.dt_s <= 0)
                throw std::runtime_error("PKD2: dt_s <= 0");
            if (hdr.n_traj == 0)
            {
                std::fclose(f);
                return {};
            }

            fseek_or_throw(f, static_cast<long>(hdr.toc_offset), SEEK_SET, "toc");
            std::vector<pkd2::TocEntry> toc(static_cast<size_t>(hdr.n_traj));
            fread_or_throw(f, toc.data(), toc.size() * sizeof(pkd2::TocEntry), "toc bytes");

            std::vector<SeedInit> out;
            out.reserve(toc.size());

            for (const auto& te : toc)
            {
                SeedInit first{}, last{};
                decode_block_first_last(hdr, te, f, first, last);

                // Choose ACTIVE and STOP
                SeedInit active = (which_point == SeedPoint::First) ? first : last;
                const SeedInit stop = (which_point == SeedPoint::First) ? last : first;

                active.has_stop     = true;
                active.lon_stop_rad = stop.lon_rad;
                active.lat_stop_rad = stop.lat_rad;
                active.t_stop_s     = stop.t_start_s;

                // Mirror ACTIVE into lon_start/lat_start for back-compat
                active.lon_start_rad = active.lon_rad;
                active.lat_start_rad = active.lat_rad;

                out.push_back(std::move(active));
            }

            std::fclose(f);
            return out;
        }
        catch (...)
        {
            std::fclose(f);
            throw;
        }
    }

    inline std::vector<SeedInit> load_seeds_from_csv_basic(const std::string& csv_path)
    {
        std::ifstream is(csv_path);
        if (!is) throw std::runtime_error("Seeds CSV: cannot open " + csv_path);

        std::string header;
        if (!std::getline(is, header))
            throw std::runtime_error("Seeds CSV: empty file " + csv_path);

        auto trim = [](std::string s){
            while (!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
            while (!s.empty() && std::isspace((unsigned char)s.back()))  s.pop_back();
            return s;
        };

        auto idx_of = [&](const std::string& name) -> int {
            std::stringstream ss(header);
            std::string tok;
            int k = 0;
            while (std::getline(ss, tok, ','))
            {
                tok = trim(tok);
                if (tok == name) return k;
                ++k;
            }
            return -1;
        };

        const int i_seed = idx_of("seed_id");
        const int i_lon  = idx_of("lon_rad");
        const int i_lat  = idx_of("lat_rad");
        const int i_t    = idx_of("t_start_s");
        if (i_seed < 0 || i_lon < 0 || i_lat < 0 || i_t < 0)
            throw std::runtime_error("Seeds CSV: expected columns: seed_id,lon_rad,lat_rad,t_start_s");

        std::vector<SeedInit> out;
        std::string line;
        while (std::getline(is, line))
        {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string f;
            int col = 0;

            SeedInit s{};
            while (std::getline(ss, f, ','))
            {
                f = trim(f);
                if (col == i_seed) s.seed_id   = std::stoull(f);
                if (col == i_lon)  s.lon_rad   = wrap_lon(std::stod(f));
                if (col == i_lat)  s.lat_rad   = std::stod(f);
                if (col == i_t)    s.t_start_s = std::stoll(f);
                ++col;
            }

            // CSV: only ACTIVE, no STOP.
            s.lon_start_rad = s.lon_rad;
            s.lat_start_rad = s.lat_rad;
            s.has_stop      = false;
            out.push_back(s);
        }
        return out;
    }
} // namespace seed_v2_detail

// =======================
// Public loaders
// =======================

inline std::vector<SeedInit> load_seeds_auto(const std::string& path,
                                             SeedPoint which_point = SeedPoint::Last)
{
    using namespace seed_v2_detail;

    if (ends_with_ci(path, ".pkd2"))
        return load_seeds_from_pkd2_zlib(path, which_point);

    if (ends_with_ci(path, ".csv"))
        return load_seeds_from_csv_basic(path);

    throw std::runtime_error("Seeds: unsupported file type (expected .pkd2 or .csv): " + path);
}

inline std::vector<SeedRec> load_seed_recs_auto(const std::string& path,
                                                SeedPoint which_point = SeedPoint::Last)
{
    return load_seeds_auto(path, which_point);
}
