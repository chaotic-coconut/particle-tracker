#ifndef DATA_PREP_UTILS_HPP
#define DATA_PREP_UTILS_HPP

#include <vector>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include <oneapi/tbb/parallel_for.h>

#include "particle_tracker/detail/ncdf_utils.hpp"
#include "particle_tracker/fixed_point/fixed_point_core.hpp"
#include "particle_tracker/interpolation_utils.hpp"

// Ensure C++17 or higher for filesystem support
#if __cplusplus < 201703L
#error "C++17 or higher is required"
#endif

// ---------------------------------------------------------------------------
// Optional debug logging (compile-time gated)
// ---------------------------------------------------------------------------
// Define DATA_PREP_DEBUG_TIME_STEPS to print time-step diagnostics.
// #define DATA_PREP_DEBUG_TIME_STEPS

// Grid struct templated on DataType
template<typename DataType>
struct Grid
{
    std::vector<DataType> lon_rad;
    std::vector<DataType> lat_rad;
    std::vector<bool> in_the_sea;
    std::size_t num_points;

    Grid() : num_points(0) {}
};

// TimeStepData struct templated on DataType and TimeType
template<typename DataType,typename TimeType=double>
struct TimeStepData
{
    TimeType time_value;
    std::map<std::string,std::vector<DataType>> variable_data;
};

// ---------------------------------------------------------------------------
// Function to initialize the grid
//
// variable_names is expected to contain:
//   [0] lon variable name
//   [1] lat variable name
//   [2] a "mask variable" used to detect sea points (NaN = land)
// ---------------------------------------------------------------------------
/**
 * @brief Initialize a lon/lat grid and a sea-mask from a netCDF file.
 *
 * The netCDF variables lon and lat are expected to be 1D arrays.
 * The mask variable is expected to match the flattened (lon x lat) grid size.
 *
 * @tparam DataType numeric type for stored grid coordinates.
 * @param file_path path to the netCDF file.
 * @param variable_names vector containing lon_var, lat_var, mask_var.
 * @return Grid<DataType> populated grid (lon_rad, lat_rad, in_the_sea).
 *
 * @throws std::runtime_error on any inconsistency or missing data.
 */
template<typename DataType>
Grid<DataType> initializeGrid(const std::string& file_path,const std::vector<std::string>& variable_names)
{
    if (variable_names.size()<3)
    {
        throw std::runtime_error("initializeGrid(): variable_names must contain lon_var, lat_var, mask_var.");
    }

    const std::string& lon_var =variable_names[0];
    const std::string& lat_var =variable_names[1];
    const std::string& mask_var=variable_names[2];

    Grid<DataType> grid;

    // Create the processor for the grid initialization
    NetCDFProcessor<DataType> processor(file_path,variable_names);

    // Retrieve variable shapes
    const auto& lon_shape=processor.getVariableShape(lon_var);
    const auto& lat_shape=processor.getVariableShape(lat_var);

    // Retrieve variable data (lon/lat are typically non-time-dependent; stored at TimeType{} key)
    const auto& lon_data=processor.getVariableData(lon_var).at(TimeType{});
    const auto& lat_data=processor.getVariableData(lat_var).at(TimeType{});

    // Use first available time slice for mask_var (works for both time-dependent and not)
    const auto& mask_map=processor.getVariableData(mask_var);
    if (mask_map.empty())
    {
        throw std::runtime_error("initializeGrid(): mask variable '"+mask_var+"' has no data.");
    }
    const auto& mask_data=mask_map.begin()->second;

    // Check if 'lon' and 'lat' are 1D arrays
    if (!(lon_shape.size()==1 && lat_shape.size()==1))
    {
        throw std::runtime_error("initializeGrid(): unsupported dimensions for '"+lon_var+"' and '"+lat_var+"'. Expected 1D lon/lat arrays.");
    }

    // Sizes of 'lon' and 'lat'
    std::size_t n_lon=lon_shape[0];
    std::size_t n_lat=lat_shape[0];

    grid.num_points=n_lon*n_lat;

    grid.lon_rad.resize(grid.num_points);
    grid.lat_rad.resize(grid.num_points);
    grid.in_the_sea.resize(grid.num_points);

    // Convert 'lon' and 'lat' to radians
    const DataType grad_to_rad_factor=DataType(pi/180.);
    std::vector<DataType> lon_rad_1d(n_lon);
    std::vector<DataType> lat_rad_1d(n_lat);

    for (std::size_t i=0;i<n_lon;i++)
    {
        lon_rad_1d[i]=lon_data[i]*grad_to_rad_factor;
    }

    for (std::size_t j=0;j<n_lat;j++)
    {
        lat_rad_1d[j]=lat_data[j]*grad_to_rad_factor;
    }

    // Create grid points by combining 'lon' and 'lat'
    for (std::size_t j=0;j<n_lat;j++)
    {
        for (std::size_t i=0;i<n_lon;i++)
        {
            std::size_t idx=j*n_lon+i;
            grid.lon_rad[idx]=lon_rad_1d[i];
            grid.lat_rad[idx]=lat_rad_1d[j];
        }
    }

    // Ensure that mask_data size matches grid.num_points
    if (mask_data.size()!=grid.num_points)
    {
        throw std::runtime_error("initializeGrid(): mask variable '"+mask_var+"' size does not match grid size.");
    }

    // Determine in_the_sea vector using mask variable (NaN = land)
    for (std::size_t idx=0;idx<grid.num_points;idx++)
    {
        const DataType m_value=mask_data[idx];
        grid.in_the_sea[idx]=!std::isnan(m_value);
    }

    return grid;
}

// ---------------------------------------------------------------------------
// Function to populate data from multiple files
// ---------------------------------------------------------------------------
/**
 * @brief Populate time-step data by reading variables from multiple files.
 *
 * This function fills all_time_steps_data with one entry per time_value.
 * For points outside the sea-mask, it writes NaN.
 *
 * @tparam DataType numeric type.
 * @tparam TimeType time type.
 * @param file_paths list of netCDF files.
 * @param grid initialized by initializeGrid().
 * @param all_time_steps_data output vector of time-step payloads.
 * @param data_variable_names list of variables to load.
 * @param check_land_points unused flag (kept for compatibility).
 */
template<typename DataType,typename TimeType=double>
void populateData(const std::vector<std::string>& file_paths,const Grid<DataType>& grid,std::vector<TimeStepData<DataType,TimeType>>& all_time_steps_data,const std::vector<std::string>& data_variable_names,bool check_land_points=false)
{
    (void)check_land_points;

    for (const auto& file_path : file_paths)
    {
        try
        {
            NetCDFProcessor<DataType,TimeType> processor(file_path,data_variable_names);

            // Processing function for data population
            auto processFunc=[&](const std::string& var_name,TimeType time_value,const std::vector<DataType>& data_vector)
            {
                // Ensure data size matches grid size
                if (data_vector.size()!=grid.num_points)
                {
                    throw std::runtime_error("Data size ("+std::to_string(data_vector.size())+
                    ") does not match grid size ("+std::to_string(grid.num_points)+").");
                }

                // Find or create TimeStepData for this time_value
                auto it=std::find_if(all_time_steps_data.begin(),all_time_steps_data.end(),
                [&](const TimeStepData<DataType,TimeType>& tsd) {return tsd.time_value==time_value;});

                if (it==all_time_steps_data.end())
                {
                    // Create a new TimeStepData
                    TimeStepData<DataType,TimeType> tsd;
                    tsd.time_value=time_value;

                    // Initialize data vectors for each variable
                    for (const auto& var : data_variable_names)
                    {
                        tsd.variable_data[var].resize(grid.num_points,std::numeric_limits<DataType>::quiet_NaN());
                    }

                    all_time_steps_data.push_back(std::move(tsd));
                    it=std::prev(all_time_steps_data.end());
                }

                // Populate data for the current variable
                auto& var_data=it->variable_data[var_name];

                for (std::size_t i=0;i<grid.num_points;i++)
                {
                    if (grid.in_the_sea[i])
                    {
                        var_data[i]=data_vector[i];
                    }
                    else
                    {
                        var_data[i]=std::numeric_limits<DataType>::quiet_NaN();
                    }
                }
            };

            // Process variables and populate data
            processor.processVariables(processFunc);
        }
        catch (const std::exception& e)
        {
            std::cerr<<"Error processing data from file '"<<file_path<<"': "<<e.what()<<std::endl;
            continue;
        }
    }
}

// Function to extract sea point indices
template<typename DataType>
std::vector<std::size_t> getSeaPointIndices(const Grid<DataType>& grid)
{
    std::vector<std::size_t> sea_point_indices;
    sea_point_indices.reserve(grid.num_points/2);

    for (std::size_t idx=0;idx<grid.num_points;idx++)
    {
        if (grid.in_the_sea[idx])
        {
            sea_point_indices.push_back(idx);
        }
    }
    return sea_point_indices;
}

template<typename TimeType>
std::vector<TimeType> generateExpectedTimeSteps(TimeType start_time,TimeType end_time,TimeType time_step_interval)
{
    std::vector<TimeType> expected_time_steps;
    for (TimeType t=start_time;t<=end_time;t+=time_step_interval)
    {
        expected_time_steps.push_back(t);
    }
    return expected_time_steps;
}

/**
 * @brief Collect time series for sea points across a set of files.
 *
 * Builds variable_time_series[var][sea_point_i][t_index].
 * The time axis is generated from earliest_time..latest_time with the provided interval.
 *
 * @tparam DataType numeric type.
 * @tparam TimeType time type (should be integer-like seconds in the pipeline).
 *
 * @throws std::runtime_error if no time steps are found.
 */
template<typename DataType,typename TimeType>
void collectTimeSeriesData(
    const Grid<DataType>& grid,
    const std::vector<std::size_t>& sea_point_indices,
    const std::vector<std::string>& variable_names,
    const std::vector<std::string>& file_paths,
    std::map<std::string, std::vector<std::vector<DataType>>>& variable_time_series,
    std::vector<TimeType>& expected_time_steps,
    TimeType time_step_interval)
{
    (void)grid;

    // Determine the start and end times
    TimeType earliest_time=std::numeric_limits<TimeType>::max();
    TimeType latest_time=std::numeric_limits<TimeType>::lowest();
    std::map<TimeType,std::string> time_to_file_map;

    // Collect all available time steps from the files
    std::set<TimeType> actual_time_steps_set;

    for (const auto& file_path : file_paths)
    {
        try
        {
            NetCDFProcessor<DataType,TimeType> processor(file_path,variable_names);

            // Get time values
            const auto& time_values=processor.getTimeValues();

            for (const auto& time_step : time_values)
            {
                actual_time_steps_set.insert(time_step);
                time_to_file_map[time_step]=file_path;
                earliest_time=std::min(earliest_time,time_step);
                latest_time=std::max(latest_time,time_step);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr<<"Error processing data from file '"<<file_path<<"': "<<e.what()<<std::endl;
            continue;
        }
    }

    if (actual_time_steps_set.empty())
    {
        throw std::runtime_error("No time steps found in data files.");
    }

#ifdef DATA_PREP_DEBUG_TIME_STEPS
    std::cerr<<"\n=== Actual time steps found in data files ===\n";
    std::cerr<<"Earliest: "<<earliest_time<<"\n";
    std::cerr<<"Latest:   "<<latest_time  <<"\n";
    std::cerr<<"All actual time steps:\n";
    for (auto t : actual_time_steps_set) std::cerr<<t<<" ";
    std::cerr<<"\n\n";
#endif

    // Generate expected time steps
    expected_time_steps=generateExpectedTimeSteps(earliest_time,latest_time,time_step_interval);
    std::size_t num_time_steps=expected_time_steps.size();

#ifdef DATA_PREP_DEBUG_TIME_STEPS
    std::cerr<<"=== Expected time steps ("<<num_time_steps<<") ===\n";
#endif

    // Initialize time series with NaNs
    for (const auto& var_name : variable_names)
    {
        variable_time_series[var_name].resize(sea_point_indices.size());
        for (auto& ts_data : variable_time_series[var_name])
        {
            ts_data.resize(num_time_steps,std::numeric_limits<DataType>::quiet_NaN());
        }
    }

    // Map actual time steps to indices in expected time steps
    std::map<TimeType,std::size_t> time_step_indices;
    for (std::size_t idx=0;idx<expected_time_steps.size();idx++)
    {
        time_step_indices[expected_time_steps[idx]]=idx;
    }

    // Load data and place it in the correct position in time series
    for (const auto& [time_step,file_path] : time_to_file_map)
    {
        try
        {
            NetCDFProcessor<DataType,TimeType> processor(file_path,variable_names);

            // Load data for this time step
            std::map<std::string,std::vector<DataType>> data_at_time;

            for (const auto& var_name : variable_names)
            {
                const auto data_vector=processor.getVariableDataAtTime(var_name,time_step);

                // Collect data for sea points
                std::vector<DataType> sea_point_data;
                sea_point_data.reserve(sea_point_indices.size());
                for (auto idx : sea_point_indices)
                {
                    sea_point_data.push_back(data_vector[idx]);
                }
                data_at_time[var_name]=std::move(sea_point_data);
            }

            // Find the index in expected_time_steps
            auto it=time_step_indices.find(time_step);
            if (it!=time_step_indices.end())
            {
                std::size_t t_idx=it->second;

                // Place data into time series
                for (const auto& var_name : variable_names)
                {
                    const auto& var_data=data_at_time.at(var_name);
                    auto& time_series=variable_time_series[var_name];

                    for (std::size_t i=0;i<sea_point_indices.size();i++)
                    {
                        time_series[i][t_idx]=var_data[i];
                    }
                }
            }
            else
            {
                std::cerr<<"Warning: Time step "<<time_step<<" not in expected time steps."<<std::endl;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr<<"Error processing data from file '"<<file_path<<"': "<<e.what()<<std::endl;
            continue;
        }
    }
}

template<typename DataType,typename TimeType=double>
void interpolateMissingData(std::vector<DataType>& time_series,const std::vector<TimeType>& time_points)
{
    std::size_t n=time_series.size();
    if (n!=time_points.size())
    {
        throw std::runtime_error("Time series and time points size mismatch.");
    }

    std::size_t start_idx=0;
    while (start_idx<n)
    {
        // Skip over existing data
        while (start_idx<n && !std::isnan(time_series[start_idx]))
        {
            start_idx++;
        }

        if (start_idx==n) break;

        // Find the end of the NaN segment
        std::size_t end_idx=start_idx;
        while (end_idx<n && std::isnan(time_series[end_idx]))
        {
            end_idx++;
        }

        // Handle edge cases
        if (start_idx==0 || end_idx==n)
        {
            // Cannot interpolate, decide how to handle
            DataType fill_value;
            if (start_idx==0 && end_idx<n)
            {
                // Missing data at the beginning, fill with first non-NaN value
                fill_value=time_series[end_idx];
            }
            else if (end_idx==n && start_idx>0)
            {
                // Missing data at the end, fill with last non-NaN value
                fill_value=time_series[start_idx-1];
            }
            else
            {
                // Entire time series is NaN
                fill_value=static_cast<DataType>(0);
            }

            for (std::size_t idx=start_idx;idx<end_idx;idx++)
            {
                time_series[idx]=fill_value;
            }
        }
        else
        {
            // Interpolate between existing data
            DataType y0=time_series[start_idx-1];
            DataType y1=time_series[end_idx];
            TimeType x0=time_points[start_idx-1];
            TimeType x1=time_points[end_idx];

            for (std::size_t idx=start_idx;idx<end_idx;idx++)
            {
                TimeType xi=time_points[idx];
                DataType yi=y0+((y1-y0)*(xi-x0)/(x1-x0));
                time_series[idx]=yi;
            }
        }

        start_idx=end_idx;
    }
}

template<typename DataType,typename TimeType=double>
void createSplinesForSeaPoints(
    std::map<std::string,std::vector<std::vector<DataType>>>&& variable_time_series,
    std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
    TimeType start_time,TimeType time_step)
{
    for (auto& [var_name,time_series_data] : variable_time_series)
    {
        std::vector<TimeSpline<DataType,TimeType>>& splines=variable_splines[var_name];
        splines.resize(time_series_data.size());

        // each i touches a unique element -> safe to move in parallel
        oneapi::tbb::parallel_for(std::size_t(0),time_series_data.size(),[&](std::size_t i)
        {
            auto& data=time_series_data[i];
            TimeSpline<DataType,TimeType> spline;
            spline.values=std::move(data);    // <-- move, no copy
            spline.createSpline(start_time,time_step);
            splines[i]=std::move(spline);
        });

        // optional: release the now-empty storage for this variable
        time_series_data.clear();
        time_series_data.shrink_to_fit();
    }

    // optional: free the whole (now moved-from) map early
    variable_time_series.clear();
    variable_time_series = {};
}

enum class EndReason : uint8_t{none=0,landed=1,neverleft=2,lifetime=3};

// --------------------------------------------------------------------
// Structure to store final trajectory results for a particle
// --------------------------------------------------------------------
template<typename DataType,typename TimeType>
struct Particle
{
    DataType init_lon;   // initial longitude (radians)
    DataType init_lat;   // initial latitude (radians)
    TimeType init_time;  // initial time (seconds since 01.01.2000 0:00)
    DataType lon;        // final longitude (radians)
    DataType lat;        // final latitude (radians)
    TimeType time;       // final time (seconds since 01.01.2000 0:00)
    EndReason reason;
};

#endif // DATA_PREP_UTILS_HPP

