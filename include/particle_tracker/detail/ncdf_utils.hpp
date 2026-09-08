#pragma once
/**
 * @file ncdf_utils.hpp
 * @brief Header-only NetCDF-C++4 helpers for loading gridded variables (optionally time-dependent).
 *
 * This header contains:
 * - particle_tracker::NetCDFDataLoader : loads one or more NetCDF files and stores variables by time.
 * - particle_tracker::NetCDFProcessor  : convenience wrapper to iterate over loaded data.
 * - particle_tracker::ncdf_utils::checkFileHeaders : fast header validation without loading full arrays.
 *
 * ## Design goals
 * - Practical: load variables into std::vector for later interpolation / tracking.
 * - Robust failure modes: prefer explicit warnings/errors over silent wrong behavior.
 * - Header-only.
 *
 * ## Important unit policy (do NOT "fix" units here)
 * - This loader does NOT normalize or convert lon/lat units.
 * - It only reads raw arrays and applies scale_factor/add_offset and missing-value handling.
 * - Unit conversion / longitude wrapping should happen in higher-level code where the expected
 *   convention is known.
 *
 * ## Missing-value policy
 * - If _FillValue or missing_value attributes exist, those are treated as missing.
 * - Otherwise (fallback) a sentinel of -30000 is treated as missing for non-coordinate variables.
 * - Missing samples are stored as quiet_NaN() in the output DataType.
 *
 * ## Time policy
 * - Reads variable "time" if present.
 * - Attempts to interpret "time:units" to convert to seconds (hours/days/seconds).
 * - Quantizes times using llround() so that common products have stable integer-second steps.
 *
 * ## Optional warning / strict mode
 * - Define PARTICLE_TRACKER_NCDF_WARN=0   to disable warnings to std::cerr.
 * - Define PARTICLE_TRACKER_NCDF_STRICT=1 to throw on certain structural problems (e.g. missing "time"
 *   while loading time-dependent variables). Default is non-strict (warn + skip).
 */

#include <algorithm>    // std::transform, std::find, std::find_if
#include <atomic>       // std::atomic_flag
#include <cctype>       // std::tolower
#include <cstddef>      // std::size_t
#include <cmath>        // std::llround
#include <exception>
#include <functional>
#include <iostream>     // std::cerr
#include <limits>       // std::numeric_limits
#include <map>
#include <set>
#include <stdexcept>    // std::runtime_error
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <netcdf>

//------------------------------------------------------------------------------
// Compile-time configuration
//------------------------------------------------------------------------------

#ifndef PARTICLE_TRACKER_NCDF_WARN
#define PARTICLE_TRACKER_NCDF_WARN 1
#endif

#ifndef PARTICLE_TRACKER_NCDF_STRICT
#define PARTICLE_TRACKER_NCDF_STRICT 0
#endif

namespace particle_tracker
{

//------------------------------------------------------------------------------
// Internal warning helpers (cheap, throttled)
//------------------------------------------------------------------------------

namespace detail
{
#if PARTICLE_TRACKER_NCDF_WARN
    inline void warn_once(const std::string& msg)
    {
        static std::atomic_flag warned = ATOMIC_FLAG_INIT;
        if (!warned.test_and_set(std::memory_order_relaxed))
        {
            std::cerr << msg << '\n';
        }
    }

    inline void warn_once_per_file(const std::string& file, const std::string& msg)
    {
        // Kept intentionally cheap: once per process. Include file name in message.
        warn_once("particle_tracker::ncdf_utils: " + file + ": " + msg);
    }
#else
    inline void warn_once(const std::string&) {}
    inline void warn_once_per_file(const std::string&, const std::string&) {}
#endif

    inline std::string to_lower(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        return s;
    }
} // namespace detail

//------------------------------------------------------------------------------
// NetCDFDataLoader
//------------------------------------------------------------------------------

/**
 * @brief Load one or more NetCDF files and store selected variables keyed by time.
 *
 * @tparam DataType numeric output type (float/double recommended)
 * @tparam TimeType time type (double by default; often integer seconds are stored after quantization)
 *
 * Data are stored as:
 *   data_[var_name][time_value] -> std::vector<DataType>
 *
 * Each vector contains the flattened variable slice for that time.
 *
 * ### Error policy
 * - Per-variable failures do not abort loadData(): they warn (optional) and continue.
 * - Accessors throw if requested data were not loaded.
 * - In strict mode (PARTICLE_TRACKER_NCDF_STRICT=1), certain structural issues throw.
 */
template<typename DataType, typename TimeType = double>
class NetCDFDataLoader
{
    public:
    using DataMap = std::map<TimeType, std::vector<DataType>>;

    NetCDFDataLoader(const std::vector<std::string>& file_paths,
                     const std::vector<std::string>& variable_names)
    : file_paths_(file_paths), variable_names_(variable_names)
    {}

    void loadData()
    {
        for (const auto& file_path : file_paths_)
        {
            loadFile(file_path);
        }
    }

    const DataMap& getVariableData(const std::string& variable_name) const
    {
        auto it = data_.find(variable_name);
        if (it != data_.end())
            return it->second;

        throw std::runtime_error("NetCDFDataLoader: variable '" + variable_name + "' not loaded.");
    }

    const std::vector<std::size_t>& getVariableShape(const std::string& variable_name) const
    {
        auto it = shapes_.find(variable_name);
        if (it != shapes_.end())
            return it->second;

        throw std::runtime_error("NetCDFDataLoader: shape for variable '" + variable_name + "' not available.");
    }

    const std::string& getVariableType(const std::string& variable_name) const
    {
        auto it = types_.find(variable_name);
        if (it != types_.end())
            return it->second;

        throw std::runtime_error("NetCDFDataLoader: type for variable '" + variable_name + "' not available.");
    }

    const std::vector<TimeType>& getTimeValues(const std::string& file_path) const
    {
        auto it = time_values_.find(file_path);
        if (it != time_values_.end())
            return it->second;

        throw std::runtime_error("NetCDFDataLoader: no time values found for file: " + file_path);
    }

    std::vector<DataType> getVariableDataAtIndex(const std::string& var_name, std::size_t index) const
    {
        auto var_it = data_.find(var_name);
        if (var_it == data_.end())
            throw std::runtime_error("NetCDFDataLoader: variable not found: " + var_name);

        if (index >= var_it->second.size())
            throw std::runtime_error("NetCDFDataLoader: index out of range for variable: " + var_name);

        auto map_iter = var_it->second.begin();
        std::advance(map_iter, static_cast<long>(index));
        return map_iter->second;
    }

    /**
     * @brief Return a netCDF variable handle (for metadata/attributes).
     *
     * Note:
     * - This opens the first file from file_paths_ and returns a variable from it.
     * - If per-file metadata are needed, open that file explicitly in caller code.
     */
    netCDF::NcVar getNetCDFVariable(const std::string& variable_name) const
    {
        if (file_paths_.empty())
            throw std::runtime_error("NetCDFDataLoader: no files configured.");

        netCDF::NcFile dataFile(file_paths_.front(), netCDF::NcFile::read);
        netCDF::NcVar var = dataFile.getVar(variable_name);
        if (var.isNull())
            throw std::runtime_error("NetCDFDataLoader: variable not found in first file: " + variable_name);

        return var;
    }

    private:
    void loadFile(const std::string& file_path)
    {
        netCDF::NcFile dataFile(file_path, netCDF::NcFile::read);

        // ------------------------------------------------------------
        // Read "time" variable (if present) and convert to seconds
        // ------------------------------------------------------------
        std::vector<TimeType> times;

        try
        {
            netCDF::NcVar time_var = dataFile.getVar("time");
            if (!time_var.isNull())
            {
                std::size_t time_size = time_var.getDim(0).getSize();
                times.resize(time_size);
                time_var.getVar(times.data());

                std::string units;
                try
                {
                    netCDF::NcVarAtt ua = time_var.getAtt("units");
                    if (!ua.isNull())
                    {
                        units.resize(256, '\0');
                        ua.getValues(units.data());
                        units.erase(std::find(units.begin(), units.end(), '\0'), units.end());
                    }
                }
                catch (...) {}

                const std::string ul = detail::to_lower(units);

                double mul = 1.0;
                if (ul.find("hour") != std::string::npos)      mul = 3600.0;
                else if (ul.find("day") != std::string::npos)  mul = 86400.0;
                else if (ul.find("sec") != std::string::npos)  mul = 1.0;
                else
                {
                    // Unknown units: keep raw but warn (once).
                    if (!units.empty())
                        detail::warn_once_per_file(file_path, "unknown time units '" + units + "'. Time values are not rescaled.");
                }

                for (auto& t : times)
                    t = static_cast<TimeType>(std::llround(static_cast<double>(t) * mul));

                time_values_[file_path] = times;
                all_time_values_.insert(times.begin(), times.end());
            }
            else
            {
#if PARTICLE_TRACKER_NCDF_STRICT
                throw std::runtime_error("missing variable 'time'");
#else
                detail::warn_once_per_file(file_path, "missing variable 'time'. Time-dependent variables may be skipped.");
#endif
            }
        }
        catch (const std::exception& e)
        {
#if PARTICLE_TRACKER_NCDF_STRICT
            throw;
#else
            detail::warn_once_per_file(file_path, std::string("error reading 'time': ") + e.what());
#endif
        }

        // ------------------------------------------------------------
        // Load requested variables
        // ------------------------------------------------------------
        for (const auto& var_name : variable_names_)
        {
            try
            {
                netCDF::NcVar var = dataFile.getVar(var_name);
                if (var.isNull())
                    throw std::runtime_error("variable not found");

                std::vector<std::size_t> shape;
                for (int i = 0; i < var.getDimCount(); i++)
                {
                    std::size_t dim_size = var.getDim(i).getSize();
                    shape.push_back(dim_size);
                }

                double scale_factor = 1.0;
                double add_offset   = 0.0;

                try
                {
                    netCDF::NcVarAtt a = var.getAtt("scale_factor");
                    if (!a.isNull()) a.getValues(&scale_factor);
                }
                catch (...) { scale_factor = 1.0; }

                try
                {
                    netCDF::NcVarAtt a = var.getAtt("add_offset");
                    if (!a.isNull()) a.getValues(&add_offset);
                }
                catch (...) { add_offset = 0.0; }

                netCDF::NcType varType = var.getType();
                const std::string typeName = varType.getName();

                if (!canConvert(varType))
                {
                    detail::warn_once_per_file(file_path,
                        "variable '" + var_name + "' has netCDF type '" + typeName + "' which cannot be converted to the template DataType.");
                    continue;
                }

                readVariableData(var, var_name, shape, scale_factor, add_offset, file_path);

                shapes_[var_name] = shape;
                types_[var_name]  = typeName;
            }
            catch (const std::exception& e)
            {
                detail::warn_once_per_file(file_path,
                    "error loading variable '" + var_name + "': " + e.what());
            }
        }
    }

    void readVariableData(netCDF::NcVar& var,
                          const std::string& var_name,
                          const std::vector<std::size_t>& shape,
                          double scale_factor,
                          double add_offset,
                          const std::string& file_path)
    {
        netCDF::NcType varType = var.getType();
        auto typeId = varType.getId();

        switch (typeId)
        {
            case NC_BYTE:   readAndTransformData<signed char>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            case NC_CHAR:   readAndTransformData<char>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            case NC_SHORT:  readAndTransformData<short>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            case NC_INT:    readAndTransformData<int>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            case NC_INT64:  readAndTransformData<long long>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            case NC_FLOAT:  readAndTransformData<float>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            case NC_DOUBLE: readAndTransformData<double>(var, var_name, shape, scale_factor, add_offset, file_path); break;
            default:
                throw std::runtime_error("unsupported netCDF variable type id");
        }
    }

    template<typename VarDataType>
    void readAndTransformData(netCDF::NcVar& var,
                              const std::string& var_name,
                              const std::vector<std::size_t>& shape,
                              double scale_factor,
                              double add_offset,
                              const std::string& file_path)
    {
        // ------------------------------------------------------------
        // Determine missing value (attribute or fallback)
        // ------------------------------------------------------------
        VarDataType missing_value = VarDataType{};
        bool has_missing_value = false;

        try
        {
            netCDF::NcVarAtt a = var.getAtt("_FillValue");
            if (!a.isNull())
            {
                a.getValues(&missing_value);
                has_missing_value = true;
            }
            else
            {
                a = var.getAtt("missing_value");
                if (!a.isNull())
                {
                    a.getValues(&missing_value);
                    has_missing_value = true;
                }
            }
        }
        catch (...) {}

        // Fallback sentinel for "data variables" (skip obvious coordinate vars)
        if (!has_missing_value && var_name != "lon" && var_name != "lat" && var_name != "depth")
        {
            missing_value = static_cast<VarDataType>(-30000);
            has_missing_value = true;
        }

        // ------------------------------------------------------------
        // Identify time dimension (if any)
        // ------------------------------------------------------------
        int time_dim_index = -1;
        for (int i = 0; i < var.getDimCount(); i++)
        {
            if (var.getDim(i).getName() == "time")
            {
                time_dim_index = i;
                break;
            }
        }

        // ------------------------------------------------------------
        // Read time-dependent variable slice-by-slice
        // ------------------------------------------------------------
        if (time_dim_index != -1)
        {
            auto it = time_values_.find(file_path);
            if (it == time_values_.end())
            {
#if PARTICLE_TRACKER_NCDF_STRICT
                throw std::runtime_error("variable '" + var_name + "' has a time dimension but no 'time' values were loaded");
#else
                detail::warn_once_per_file(file_path,
                    "variable '" + var_name + "' has a time dimension but no 'time' values were loaded. Skipping.");
                return;
#endif
            }

            const auto& times = it->second;
            const std::size_t time_size = shape[static_cast<std::size_t>(time_dim_index)];

            if (times.size() < time_size)
            {
#if PARTICLE_TRACKER_NCDF_STRICT
                throw std::runtime_error("time vector shorter than time dimension for variable '" + var_name + "'");
#else
                detail::warn_once_per_file(file_path,
                    "time vector shorter than time dimension for variable '" + var_name + "'. Skipping.");
                return;
#endif
            }

            std::vector<std::size_t> start(static_cast<std::size_t>(var.getDimCount()), 0);
            std::vector<std::size_t> count = shape;
            count[static_cast<std::size_t>(time_dim_index)] = 1;

            for (std::size_t t = 0; t < time_size; t++)
            {
                start[static_cast<std::size_t>(time_dim_index)] = t;
                TimeType time_value = times[t];

                std::size_t data_size = 1;
                for (std::size_t dim_size : count)
                    data_size *= dim_size;

                std::vector<VarDataType> raw_data(data_size);
                var.getVar(start, count, raw_data.data());

                std::vector<DataType> transformed_data(data_size);

                for (std::size_t i = 0; i < data_size; i++)
                {
                    if (has_missing_value && raw_data[i] == missing_value)
                    {
                        transformed_data[i] = std::numeric_limits<DataType>::quiet_NaN();
                    }
                    else
                    {
                        transformed_data[i] =
                            static_cast<DataType>(raw_data[i]) * static_cast<DataType>(scale_factor) +
                            static_cast<DataType>(add_offset);
                    }
                }

                data_[var_name][time_value] = std::move(transformed_data);
            }
        }
        else
        {
            // ------------------------------------------------------------
            // Read variable without time dimension
            // ------------------------------------------------------------
            std::size_t total_size = 1;
            for (std::size_t dim_size : shape)
                total_size *= dim_size;

            std::vector<VarDataType> raw_data(total_size);
            var.getVar(raw_data.data());

            std::vector<DataType> transformed_data(total_size);

            for (std::size_t i = 0; i < total_size; i++)
            {
                if (has_missing_value && raw_data[i] == missing_value)
                {
                    transformed_data[i] = std::numeric_limits<DataType>::quiet_NaN();
                }
                else
                {
                    transformed_data[i] =
                        static_cast<DataType>(raw_data[i]) * static_cast<DataType>(scale_factor) +
                        static_cast<DataType>(add_offset);
                }
            }

            TimeType time_value = TimeType{};
            data_[var_name][time_value] = std::move(transformed_data);
        }
    }

    bool canConvert(const netCDF::NcType& varType) const
    {
        return std::is_arithmetic<DataType>::value && (
            varType == netCDF::ncByte   ||
            varType == netCDF::ncChar   ||
            varType == netCDF::ncShort  ||
            varType == netCDF::ncInt    ||
            varType == netCDF::ncInt64  ||
            varType == netCDF::ncFloat  ||
            varType == netCDF::ncDouble
        );
    }

    std::vector<std::string> file_paths_;
    std::vector<std::string> variable_names_;

    std::map<std::string, std::map<TimeType, std::vector<DataType>>> data_;
    std::map<std::string, std::vector<std::size_t>> shapes_;
    std::map<std::string, std::string> types_;
    std::map<std::string, std::vector<TimeType>> time_values_;
    std::set<TimeType> all_time_values_;
};

//------------------------------------------------------------------------------
// NetCDFProcessor
//------------------------------------------------------------------------------

/**
 * @brief Convenience wrapper to load a single file and iterate over selected variables and times.
 */
template<typename DataType, typename TimeType = double>
class NetCDFProcessor
{
    public:
    using DataMap = std::map<TimeType, std::vector<DataType>>;
    using VariableProcessingFunc =
        std::function<void(const std::string&, TimeType, const std::vector<DataType>&)>;

    NetCDFProcessor(const std::string& file_path,
                    const std::vector<std::string>& variable_names)
    : file_path_(file_path),
      variable_names_(variable_names),
      data_loader_({file_path}, variable_names)
    {
        try
        {
            data_loader_.loadData();
        }
        catch (const std::exception& e)
        {
            std::cerr << "NetCDFProcessor: error loading file '" << file_path_ << "': " << e.what() << '\n';
            throw;
        }
    }

    void processVariables(const VariableProcessingFunc& processFunc)
    {
        try
        {
            for (const auto& var_name : variable_names_)
            {
                const auto& var_data = data_loader_.getVariableData(var_name);
                for (const auto& kv : var_data)
                {
                    processFunc(var_name, kv.first, kv.second);
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "NetCDFProcessor: error processing file '" << file_path_ << "': " << e.what() << '\n';
            throw;
        }
    }

    netCDF::NcVar getNetCDFVariable(const std::string& variable_name) const
    {
        return data_loader_.getNetCDFVariable(variable_name);
    }

    const std::vector<std::size_t>& getVariableShape(const std::string& variable_name) const
    {
        return data_loader_.getVariableShape(variable_name);
    }

    const std::string& getVariableType(const std::string& variable_name) const
    {
        return data_loader_.getVariableType(variable_name);
    }

    const std::vector<TimeType>& getTimeValues() const
    {
        return data_loader_.getTimeValues(file_path_);
    }

    const typename NetCDFDataLoader<DataType, TimeType>::DataMap&
    getVariableData(const std::string& variable_name) const
    {
        return data_loader_.getVariableData(variable_name);
    }

    /*std::vector<DataType> getVariableDataAtTime(const std::string& var_name, TimeType time_value)
    {
        std::size_t time_index = getTimeIndex(time_value);
        return data_loader_.getVariableDataAtIndex(var_name, time_index);
    }*/

    std::vector<DataType> getVariableDataAtTime(const std::string& var_name, TimeType time_value)
    {
        const auto& var_data = data_loader_.getVariableData(var_name);
    
        const TimeType key = static_cast<TimeType>(std::llround(static_cast<double>(time_value)));
    
        auto it = var_data.find(key);
        if (it == var_data.end())
            throw std::runtime_error("NetCDFProcessor: requested time slice not loaded for variable '" +
                                     var_name + "' at time=" + std::to_string(static_cast<long long>(key)));
    
        return it->second; // returns a copy (your current API returns by value anyway)
    }

    const std::vector<DataType>& getVariableDataAtTimeRef(const std::string& var_name, TimeType time_value) const
    {
        const auto& var_data = data_loader_.getVariableData(var_name);
        const TimeType key = static_cast<TimeType>(std::llround(static_cast<double>(time_value)));
    
        auto it = var_data.find(key);
        if (it == var_data.end())
            throw std::runtime_error("NetCDFProcessor: requested time slice not loaded...");
    
        return it->second;
    }

    private:
    std::string file_path_;
    std::vector<std::string> variable_names_;
    NetCDFDataLoader<DataType, TimeType> data_loader_;

    std::size_t getTimeIndex(TimeType time_value)
    {
        const auto& time_values = data_loader_.getTimeValues(file_path_);
        const long long q = std::llround(static_cast<double>(time_value));

        auto it = std::find_if(time_values.begin(), time_values.end(),
                               [q](TimeType t){ return std::llround(static_cast<double>(t)) == q; });

        if (it != time_values.end())
            return static_cast<std::size_t>(std::distance(time_values.begin(), it));

        throw std::runtime_error("NetCDFProcessor: requested time value not found in loaded data.");
    }
};

//------------------------------------------------------------------------------
// ncdf_utils helper functions
//------------------------------------------------------------------------------

namespace ncdf_utils
{
/**
 * @brief Validate netCDF file headers (metadata-only; does not read full arrays).
 */
inline std::vector<std::string> checkFileHeaders(const std::string& fileName,
                                                 const std::vector<std::string>& requiredVariables,
                                                 int expectedTimeCount)
{
    std::vector<std::string> errors;

    try
    {
        netCDF::NcFile dataFile(fileName, netCDF::NcFile::read);

        netCDF::NcVar timeVar = dataFile.getVar("time");
        if (timeVar.isNull())
        {
            errors.push_back("Missing variable 'time'");
        }
        else
        {
            if (timeVar.getDimCount() < 1)
            {
                errors.push_back("Variable 'time' has no dimensions.");
            }
            else
            {
                std::size_t timeSize = timeVar.getDim(0).getSize();
                if (timeSize != static_cast<std::size_t>(expectedTimeCount))
                {
                    errors.push_back("Variable 'time' dimension size is " +
                                     std::to_string(timeSize) + " but expected " +
                                     std::to_string(expectedTimeCount));
                }
            }
        }

        for (const auto& varName : requiredVariables)
        {
            if (varName == "time") continue;

            netCDF::NcVar var = dataFile.getVar(varName);
            if (var.isNull())
            {
                errors.push_back("Missing variable '" + varName + "'");
                continue;
            }

            int dimCount = var.getDimCount();
            if (dimCount == 0)
            {
                errors.push_back("Variable '" + varName + "' has no dimensions.");
                continue;
            }

            for (int i = 0; i < dimCount; i++)
            {
                std::size_t dSize = var.getDim(i).getSize();
                if (dSize == 0)
                {
                    errors.push_back("Variable '" + varName + "' has dimension " +
                                     std::to_string(i) + " equal to 0.");
                }
            }

            if (dimCount > 1)
            {
                std::size_t timeDim = var.getDim(0).getSize();
                if (timeDim != static_cast<std::size_t>(expectedTimeCount))
                {
                    errors.push_back("Variable '" + varName + "' has time dimension " +
                                     std::to_string(timeDim) + " but expected " +
                                     std::to_string(expectedTimeCount));
                }
            }
        }
    }
    catch (const std::exception& e)
    {
        errors.push_back("Error processing file '" + fileName + "': " + std::string(e.what()));
    }

    return errors;
}
} // namespace ncdf_utils

} // namespace particle_tracker
