/********************************************************************
 *  fixed_point_pkd2.hpp  -  per-trajectory compressed storage (PKD-2)
 *  Requires: fixed_point_core.hpp
 *  Define PK_IO_USE_ZLIB and link with zlib.
 ********************************************************************/

// NOTE ON PKD2 BLOCK SIZE LIMITS
//
// PKD2 stores each trajectory as one independently-compressed block.
// The TOC stores `raw_size` and `comp_size` as uint32_t, therefore:
//
//   * Each individual trajectory block must be < 4 GiB uncompressed.
//   * Each individual trajectory block must be < 4 GiB compressed.
//
// Additionally, the zlib single-shot APIs used here (compress2/uncompress)
// take sizes as uLong/uLongf (platform-dependent). On some platforms
// (notably Windows), this is effectively a 4 GiB limit as well.
// If you ever need trajectories larger than this, you must split them
// into multiple PKD2 blocks (i.e., change the file format), or implement
// streaming compression and store 64-bit sizes in the TOC.

#pragma once
#include "fixed_point_core.hpp" // fp_pi, quant_coord, wrap_lon, ...
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#ifdef PK_IO_USE_ZSTD
#error "PK_IO_USE_ZSTD is unsupported; PKD2 files use zlib compression"
#endif
#ifndef PK_IO_USE_ZLIB
#error "Define PK_IO_USE_ZLIB before including fixed_point_pkd2.hpp"
#else
#include <zlib.h>
#endif

inline constexpr int kDefaultLevel = 6;

namespace pkd2 {

/* --------------------------- file layout --------------------------- */
#pragma pack(push, 1)
struct Header {
  char magic[8] = {'P', 'K', 'D', '2', 0, 0, 0, 0};
  uint32_t version = 3; // bump
  uint8_t little = 1;
  uint8_t _r1{0};
  uint16_t _r2{0};
  uint32_t pos_scale_ticks_per_rad = 1'000'000;
  int64_t t0_epoch_s = 0;
  int32_t dt_s = 0; // global base step (e.g., 3h)
  uint32_t _r3{0};
  uint64_t n_traj = 0;
  uint64_t toc_offset = 0;
};

// header of the uncompressed per-trajectory block
struct BlockHeader {
  uint64_t seed_id;
  int32_t start_index; // index of first stored sample on the global 3h grid
  uint16_t stride;     // >=1  (number of base steps between stored samples)
  uint16_t _pad{0};    // keep 4-byte alignment
  uint32_t n_points;   // >=1 number of stored samples
  int32_t lon0_q;      // ticks = radians * pos_scale_ticks_per_rad
  int32_t lat0_q;      // ticks = radians * pos_scale_ticks_per_rad
  // then (n_points-1) dlon_q[], then (n_points-1) dlat_q[]  (int32 each)
};

// TOC entry at EOF
struct TocEntry {
  uint64_t seed_id;
  uint64_t file_offset;
  uint32_t comp_size;
  uint32_t raw_size;
  uint32_t n_points;
  int32_t start_index;
  uint16_t stride; // >=1
  uint16_t _pad{0};
};
#pragma pack(pop)

/* --------------------------- writer -------------------------------- */
class Writer {
public:
  // open for write; header is written immediately
  Writer(const std::string &path, int64_t t0_epoch_s, int32_t dt_s,
         uint32_t pos_scale = 1'000'000, int compression_level = kDefaultLevel)
      : lvl_(compression_level), pos_scale_(pos_scale) {
    f_ = std::fopen(path.c_str(), "wb+");
    if (!f_)
      throw std::runtime_error("pkd2::Writer: cannot open " + path);
    std::memset(&hdr_, 0, sizeof(hdr_));
    hdr_.magic[0] = 'P';
    hdr_.magic[1] = 'K';
    hdr_.magic[2] = 'D';
    hdr_.magic[3] = '2';
    hdr_.version = 3;
    hdr_.little = 1;
    hdr_.pos_scale_ticks_per_rad = pos_scale_;
    hdr_.t0_epoch_s = t0_epoch_s;
    hdr_.dt_s = dt_s;
    write_at(0, &hdr_, sizeof(hdr_));
    std::fseek(f_, 0, SEEK_END);
  }

  ~Writer() {
    try {
      close();
    } catch (...) {
    }
  }

  // Append a trajectory given lon/lat in **radians** (the native state)
  // Optional stride lets to decimate (e.g. stride=8 turns 3h to 24h).
  template <class Real>
  void add_traj_rad(uint64_t seed_id, const Real *lon_rad, const Real *lat_rad,
                    uint32_t n_points, int32_t start_index,
                    uint32_t stride = 1) {
    static_assert(std::is_floating_point_v<Real>, "Real must be float/double");
    if (n_points == 0)
      return;

    const bool need_copy = stride > 1 || !std::is_same_v<Real, double>;

    if (need_copy) {
      tmp_lon_.clear();
      tmp_lat_.clear();
      tmp_lon_.reserve((n_points + stride - 1) / stride);
      tmp_lat_.reserve(tmp_lon_.capacity());
      for (uint32_t i = 0; i < n_points; i += stride) {
        tmp_lon_.push_back(static_cast<double>(lon_rad[i]));
        tmp_lat_.push_back(static_cast<double>(lat_rad[i]));
      }
      emit_block(seed_id, tmp_lon_.data(), tmp_lat_.data(),
                 static_cast<uint32_t>(tmp_lon_.size()), start_index,
                 static_cast<uint16_t>(stride));
    } else {
      // Real == double and stride == 1
      emit_block(seed_id, reinterpret_cast<const double *>(lon_rad),
                 reinterpret_cast<const double *>(lat_rad), n_points,
                 start_index, /*stride*/ 1);
    }
  }

  void close() {
    if (!f_)
      return;
    // write TOC
    long pos = std::ftell(f_);
    if (pos < 0)
      throw std::runtime_error("pkd2::Writer: ftell failed");
    if (!toc_.empty()) {
      size_t nw = std::fwrite(toc_.data(), sizeof(TocEntry), toc_.size(), f_);
      if (nw != toc_.size())
        throw std::runtime_error("pkd2::Writer: fwrite TOC failed");
    }
    std::fflush(f_);
    hdr_.n_traj = static_cast<uint64_t>(toc_.size());
    hdr_.toc_offset = static_cast<uint64_t>(pos);
    write_at(0, &hdr_, sizeof(hdr_));
    std::fflush(f_);
    std::fclose(f_);
    f_ = nullptr;
  }

private:
  FILE *f_{nullptr};
  Header hdr_{};
  int lvl_;
  uint32_t pos_scale_;
  std::vector<TocEntry> toc_;
  std::vector<double> tmp_lon_, tmp_lat_;
  std::vector<int32_t> dlon_, dlat_;
  std::vector<unsigned char> raw_, comp_;

  static inline double clamp_lat(double r) {
    const double hp = fp_pi * .5;
    return (r < -hp) ? -hp : (r > hp ? hp : r);
  }
  static inline double unwrap_pi(double r) {
    // (-pi,pi]
    while (r <= -fp_pi)
      r += 2 * fp_pi;
    while (r > fp_pi)
      r -= 2 * fp_pi;
    return r;
  }
  static inline double unwrap_follow(double prev, double r) {
    r = unwrap_pi(r);
    double d = r - prev;
    if (d > fp_pi)
      r -= 2 * fp_pi;
    if (d < -fp_pi)
      r += 2 * fp_pi;
    return r;
  }
  inline int32_t to_q(double radians) const {
    long long v = llround(radians * static_cast<double>(pos_scale_));
    if (v > std::numeric_limits<int32_t>::max())
      v = std::numeric_limits<int32_t>::max();
    if (v < std::numeric_limits<int32_t>::min())
      v = std::numeric_limits<int32_t>::min();
    return static_cast<int32_t>(v);
  }

  void emit_block(uint64_t seed_id, const double *lon_r, const double *lat_r,
                  uint32_t n, int32_t start_index, uint16_t stride) {
    // build uncompressed block [BlockHeader | dlon[n-1] | dlat[n-1]]
    BlockHeader bh{};
    bh.seed_id = seed_id;
    bh.start_index = start_index;
    bh.stride = stride ? stride : 1;
    bh.n_points = n;

    double lon0 = unwrap_pi(lon_r[0]);
    double lat0 = clamp_lat(lat_r[0]);
    bh.lon0_q = to_q(lon0);
    bh.lat0_q = to_q(lat0);

    dlon_.assign(n > 1 ? n - 1 : 0, 0);
    dlat_.assign(n > 1 ? n - 1 : 0, 0);

    double prevL = lon0, prevB = lat0;
    for (uint32_t i = 1; i < n; i++) {
      double lr = unwrap_follow(prevL, lon_r[i]);
      double br = clamp_lat(lat_r[i]);
      dlon_[i - 1] = to_q(lr - prevL);
      dlat_[i - 1] = to_q(br - prevB);
      prevL = lr;
      prevB = br;
    }

    const size_t raw_bytes = sizeof(BlockHeader) +
                             dlon_.size() * sizeof(int32_t) +
                             dlat_.size() * sizeof(int32_t);

    // PKD2 format constraint: TOC stores sizes as uint32_t.
    if (raw_bytes > static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
      throw std::runtime_error("pkd2: raw block exceeds 4GiB (format limit)");

    raw_.resize(raw_bytes);
    std::memcpy(raw_.data(), &bh, sizeof(bh));
    if (!dlon_.empty()) {
      std::memcpy(raw_.data() + sizeof(bh), dlon_.data(),
                  dlon_.size() * sizeof(int32_t));

      std::memcpy(raw_.data() + sizeof(bh) + dlon_.size() * sizeof(int32_t),
                  dlat_.data(), dlat_.size() * sizeof(int32_t));
    }

    // compress to memory
    uint32_t comp_sz = 0;

    {
      // zlib APIs here take uLong/uLongf (platform-dependent width).
      if (raw_bytes > static_cast<size_t>(std::numeric_limits<uLong>::max()))
        throw std::runtime_error("pkd2: raw block too large for zlib uLong");

      const uLongf bound = compressBound(static_cast<uLong>(raw_bytes));

      if (static_cast<size_t>(bound) >
          static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
        throw std::runtime_error(
            "pkd2: zlib compressBound exceeds 4GiB (format limit)");

      comp_.resize(static_cast<size_t>(bound));

      uLongf z = bound;
      const int rc = compress2(comp_.data(), &z,
                               reinterpret_cast<const Bytef *>(raw_.data()),
                               static_cast<uLong>(raw_bytes), lvl_);
      if (rc != Z_OK)
        throw std::runtime_error("pkd2: compress2 failed");

      if (static_cast<size_t>(z) >
          static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
        throw std::runtime_error(
            "pkd2: compressed block exceeds 4GiB (format limit)");

      comp_sz = static_cast<uint32_t>(z);
    }

    // write to file, remember offset
    long off = std::ftell(f_);
    if (off < 0)
      throw std::runtime_error("pkd2: ftell failed");

    const size_t nw = std::fwrite(comp_.data(), 1, comp_sz, f_);
    if (nw != comp_sz)
      throw std::runtime_error("pkd2: fwrite failed");

    TocEntry te{};
    te.seed_id = seed_id;
    te.file_offset = static_cast<uint64_t>(off);
    te.comp_size = comp_sz;
    te.raw_size = static_cast<uint32_t>(raw_bytes);
    te.n_points = n;
    te.start_index = start_index;
    te.stride = bh.stride;

    toc_.push_back(te);
  }

  void write_at(long off, const void *p, size_t n) {
    if (std::fseek(f_, off, SEEK_SET) != 0)
      throw std::runtime_error("fseek failed");
    size_t nw = std::fwrite(p, 1, n, f_);
    if (nw != n)
      throw std::runtime_error("header write failed");
    std::fseek(f_, 0, SEEK_END);
  }
};

/* ---------------------------- reader ------------------------------- */
class Reader {
public:
  explicit Reader(const std::string &path) {
    f_ = std::fopen(path.c_str(), "rb");
    if (!f_)
      throw std::runtime_error("pkd2::Reader: cannot open " + path);
    read_at(0, &hdr_, sizeof(hdr_));
    if (hdr_.magic[0] != 'P' || hdr_.magic[1] != 'K' || hdr_.magic[2] != 'D' ||
        hdr_.magic[3] != '2')
      throw std::runtime_error("pkd2: bad magic");
    if (hdr_.version != 3)
      throw std::runtime_error("pkd2: version != 3");
    // read TOC
    toc_.resize(hdr_.n_traj);
    read_at(static_cast<long>(hdr_.toc_offset), toc_.data(),
            toc_.size() * sizeof(TocEntry));
    for (size_t i = 0; i < toc_.size(); ++i)
      index_[toc_[i].seed_id] = static_cast<uint32_t>(i);
  }
  ~Reader() {
    if (f_)
      std::fclose(f_);
  }

  bool has(uint64_t seed_id) const {
    return index_.find(seed_id) != index_.end();
  }

  // Read one trajectory; outputs in **radians** (to match sim)
  void read_traj(uint64_t seed_id, std::vector<double> &lon_r,
                 std::vector<double> &lat_r, std::vector<int64_t> &t_s) const {
    auto it = index_.find(seed_id);
    if (it == index_.end())
      throw std::runtime_error("pkd2: seed_id not found");
    const TocEntry &te = toc_[it->second];

    // read compressed block
    comp_.resize(te.comp_size);
    read_at(static_cast<long>(te.file_offset), comp_.data(), comp_.size());

    // decompress
    raw_.resize(te.raw_size);
    uLongf out = static_cast<uLongf>(raw_.size());
    int rc = uncompress(raw_.data(), &out, comp_.data(),
                        static_cast<uLong>(comp_.size()));
    if (rc != Z_OK || out != raw_.size())
      throw std::runtime_error("uncompress failed");
    // parse block
    if (raw_.size() < sizeof(BlockHeader))
      throw std::runtime_error("pkd2: raw too small");
    const BlockHeader *bh = reinterpret_cast<const BlockHeader *>(raw_.data());
    const uint16_t s = bh->stride ? bh->stride : 1;
    if (bh->n_points != te.n_points)
      throw std::runtime_error("pkd2: n_points mismatch");

    const int32_t *dlon =
        reinterpret_cast<const int32_t *>(raw_.data() + sizeof(BlockHeader));
    const int32_t *dlat = dlon + (bh->n_points > 0 ? bh->n_points - 1 : 0);

    lon_r.resize(bh->n_points);
    lat_r.resize(bh->n_points);
    t_s.resize(bh->n_points);

    // reconstruct
    const double inv_scale =
        1. / static_cast<double>(hdr_.pos_scale_ticks_per_rad);
    double L = static_cast<double>(bh->lon0_q) * inv_scale;
    double B = static_cast<double>(bh->lat0_q) * inv_scale;
    lon_r[0] = L;
    lat_r[0] = B;
    t_s[0] =
        hdr_.t0_epoch_s + static_cast<int64_t>(bh->start_index) * hdr_.dt_s;

    for (uint32_t i = 1; i < bh->n_points; i++) {
      L += static_cast<double>(dlon[i - 1]) * inv_scale;
      B += static_cast<double>(dlat[i - 1]) * inv_scale;
      lon_r[i] = L;
      lat_r[i] = B;
      t_s[i] = t_s[0] + static_cast<int64_t>(i) * hdr_.dt_s * s;
    }
  }

  const Header &header() const { return hdr_; }
  const std::vector<TocEntry> &toc() const { return toc_; }

private:
  FILE *f_{nullptr};
  Header hdr_{};
  std::vector<TocEntry> toc_;
  std::unordered_map<uint64_t, uint32_t> index_;
  mutable std::vector<unsigned char> comp_, raw_;

  void read_at(long off, void *p, size_t n) const {
    if (std::fseek(f_, off, SEEK_SET) != 0)
      throw std::runtime_error("fseek failed");
    size_t nr = std::fread(p, 1, n, f_);
    if (nr != n)
      throw std::runtime_error("fread failed");
  }
};

} // namespace pkd2
