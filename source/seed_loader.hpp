// seed_loader.hpp
#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <random>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include "fixed_point_core.hpp"  // for wrap_lon and pi

inline constexpr double POS_TICKS_PER_RAD=1'000'000.0; // micro-rad per tick

struct SeedInit
{
	uint64_t seed_id;
	double   lon_rad;   // start lon (Hawaii side), radians in (-pi, pi]
	double   lat_rad;   // start lat, radians
	int64_t  t_start_s; // start time (Unix seconds since 1970-01-01 00:00:00 UTC (downstream converts to t2000))
};

// Deterministic 64-bit mixer (splitmix64)
static inline uint64_t mix64(uint64_t x)
{
	x+=0x9e3779b97f4a7c15ULL;
	x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
	x=(x^(x>>27))*0x94d049bb133111ebULL;
	return x^(x>>31);
}

/* ---------------- Arrow/Parquet path ---------------- */
#ifdef USE_ARROW
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

// Read selected columns from a single Parquet file and return sampled seeds.
// - sample_frac in (0,1] keeps that fraction deterministically by hashing seed_id.
// - sample_max > 0 caps the number of returned seeds.
// - seed_offset lets ``shift'' the hash for independent sub-samples.
inline std::vector<SeedInit>
load_seeds_from_parquet(const std::string& path,
                        double sample_frac=1.,
                        uint64_t sample_max=0,
                        uint64_t seed_offset=0)
{
	if (sample_frac<=0.) return {};
	if (sample_frac>1.)  sample_frac=1.0;

	ARROW_ASSIGN_OR_RAISE(auto infile,arrow::io::ReadableFile::Open(path));
	std::unique_ptr<parquet::arrow::FileReader> reader;
	PARQUET_THROW_NOT_OK(parquet::arrow::OpenFile(infile,arrow::default_memory_pool(),&reader));

	// Only pull what we need
	std::shared_ptr<arrow::Table> table;
	// If libarrow is new enough, use ReadTable with column projection by names.
	// For broad compatibility, read full table and pick columns by name.
	PARQUET_THROW_NOT_OK(reader->ReadTable(&table));

	auto get=[&](const char* name)->std::shared_ptr<arrow::ChunkedArray>
	{
		auto col=table->GetColumnByName(name);
		if (!col) throw std::runtime_error(std::string("Parquet: missing column '")+name+"'");
		return col;
	};

	auto col_seed=get("seed_id");     // uint64
	auto col_lonq=get("lon_hi_q");    // int64 (micro-rad)
	auto col_latq=get("lat_hi_q");    // int64 (micro-rad)
	auto col_ts  =get("t_hi_s");      // int64 (simulation seconds)

	// Iterate by record batch to keep memory locality and column alignment
	arrow::TableBatchReader tbr(*table);
	tbr.set_chunksize(256'000);

	std::vector<SeedInit> out;
	std::shared_ptr<arrow::RecordBatch> batch;

	const long double TH=(long double)std::numeric_limits<uint64_t>::max()*static_cast<long double>(sample_frac);

	while (true)
	{
		ARROW_ASSIGN_OR_RAISE(batch,tbr.Next());
		if (!batch) break;

		auto a_seed=std::static_pointer_cast<arrow::UInt64Array>(batch->GetColumnByName("seed_id"));
		auto a_lonq=std::static_pointer_cast<arrow::Int64Array >(batch->GetColumnByName("lon_hi_q"));
		auto a_latq=std::static_pointer_cast<arrow::Int64Array >(batch->GetColumnByName("lat_hi_q"));
		auto a_ts  =std::static_pointer_cast<arrow::Int64Array >(batch->GetColumnByName("t_hi_s"));
		if (!a_seed || !a_lonq || !a_latq || !a_ts)
			throw std::runtime_error("Parquet: column type mismatch");

		const int64_t n=batch->num_rows();
		out.reserve(out.size()+static_cast<size_t>(n*sample_frac*1.05));

		for (int64_t i=0;i<n;++i)
		{
			if (a_seed->IsNull(i) || a_lonq->IsNull(i) || a_latq->IsNull(i) || a_ts->IsNull(i)) continue;

			uint64_t sid=a_seed->Value(i);
			uint64_t h  =mix64(sid^seed_offset);
			if ((long double)h>TH) continue; // sampled out

			double lon_rad=static_cast<double>(a_lonq->Value(i)/POS_TICKS_PER_RAD);
			double lat_rad=static_cast<double>(a_latq->Value(i)/POS_TICKS_PER_RAD);

			// Wrap lon into (-pi,pi] to match your simulation
			lon_rad=wrap_lon(lon_rad);

			out.push_back(SeedInit{
						sid,
						lon_rad,
						lat_rad,
						static_cast<int64_t>(a_ts->Value(i))
			});

			if (sample_max && out.size()>=sample_max) return out;
		}
	}
	return out;
}
#else
/* ---------------- CSV fallback ---------------- */
// If Arrow is not linked, accept a CSV with the same column names.
inline std::vector<SeedInit>
load_seeds_from_parquet(const std::string& csv_path,
                        double sample_frac=1.,
                        uint64_t sample_max=0,
                        uint64_t seed_offset=0)
{
	std::ifstream is(csv_path);
	if (!is) throw std::runtime_error("cannot open CSV: "+csv_path);
	std::string header;
	std::getline(is,header);

	// crude column index finder
	auto idx=[&](const std::string& name)
	{
		std::stringstream ss(header);
		std::string tok;int k=0;
		while (std::getline(ss,tok,','))
		{
			if (tok==name) return k;
			++k;
		}
		return -1;
	};
	int i_seed=idx("seed_id"),i_lonq=idx("lon_hi_q"),i_latq=idx("lat_hi_q"),i_ts=idx("t_hi_s");
	if (i_seed < 0 || i_lonq < 0 || i_latq < 0 || i_ts < 0)
		throw std::runtime_error("CSV: header missing required columns");

	const long double TH=(long double)std::numeric_limits<uint64_t>::max()*static_cast<long double>(sample_frac);

	std::vector<SeedInit> out;
	std::string line;
	while (std::getline(is,line))
	{
		std::stringstream ss(line);
		std::string f; int col=0;
		uint64_t sid=0;long long lonq=0,latq=0,ts=0;

		while (std::getline(ss,f,','))
		{
			if (col==i_seed) sid =std::stoull(f);
			if (col==i_lonq) lonq=std::stoll(f);
			if (col==i_latq) latq=std::stoll(f);
			if (col==i_ts)   ts  =std::stoll(f);
			++col;
		}

		uint64_t h=mix64(sid^seed_offset);
		if ((long double)h>TH) continue;

		double lon_rad=wrap_lon(lonq/POS_TICKS_PER_RAD);
		double lat_rad=latq/POS_TICKS_PER_RAD;

		out.push_back(SeedInit{sid,lon_rad,lat_rad,ts});
		if (sample_max && out.size()>=sample_max) break;
	}
	return out;
}
#endif

