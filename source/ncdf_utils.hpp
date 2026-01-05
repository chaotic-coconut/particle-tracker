#ifndef NCDF_UTILS_HPP
#define NCDF_UTILS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <exception>
#include <limits>
#include <netcdf>

template<typename DataType,typename TimeType=double>
class NetCDFDataLoader
{
	public:
	// Constructor accepts file paths and variable names to load
	NetCDFDataLoader(const std::vector<std::string>& file_paths,const std::vector<std::string>& variable_names)
        : file_paths_(file_paths),variable_names_(variable_names) {}

	using DataMap=std::map<TimeType,std::vector<DataType>>;

	// Load data from all files
	void loadData()
	{
		for (const auto& file_path : file_paths_)
		{
			loadFile(file_path);
		}
	}

	// Get variable data
	const std::map<TimeType,std::vector<DataType>>& getVariableData(const std::string& variable_name) const
	{
		auto it=data_.find(variable_name);
		if (it!=data_.end())
		{
			return it->second;
		}
		else
		{
			throw std::runtime_error("Variable '"+variable_name+"' not loaded.");
		}
	}

	// Get the shape of a variable
	const std::vector<std::size_t>& getVariableShape(const std::string& variable_name) const
	{
		auto it=shapes_.find(variable_name);
		if (it!=shapes_.end())
		{
			return it->second;
		}
		else
		{
			throw std::runtime_error("Shape for variable '"+ variable_name+"' not available.");
		}
	}

	// Get the type name of a variable
	const std::string& getVariableType(const std::string& variable_name) const
	{
		auto it=types_.find(variable_name);
		if (it!=types_.end())
		{
			return it->second;
		}
		else
		{
			throw std::runtime_error("Type for variable '"+variable_name+"' not available.");
		}
	}

	const std::vector<TimeType>& getTimeValues(const std::string& file_path) const
	{
		auto it=time_values_.find(file_path);
		if (it!=time_values_.end())
		{
			return it->second;  // Return the vector of times
		}
		else
		{
			throw std::runtime_error("No time values found for file: "+file_path);
		}
	}

	std::vector<DataType> getVariableDataAtIndex(const std::string& var_name,size_t index) const
	{
		auto var_it=data_.find(var_name);
		if (var_it==data_.end())
		{
			throw std::runtime_error("Variable not found: "+var_name);
		}

		if (index>=var_it->second.size())
		{
			throw std::runtime_error("Index out of range for variable "+var_name);
		}

		auto map_iter=var_it->second.begin();
		std::advance(map_iter,index);

		return map_iter->second;
	}

	private:
	void loadFile(const std::string& file_path)
	{
		netCDF::NcFile dataFile(file_path,netCDF::NcFile::read);

		// Read the 'time' variable
		std::vector<TimeType> times;
		try
		{
			netCDF::NcVar time_var=dataFile.getVar("time");
			if (!time_var.isNull())
			{
				std::size_t time_size=time_var.getDim(0).getSize();
				times.resize(time_size);
				time_var.getVar(times.data());

				// Convert to seconds based on units
				std::string units;
				try
				{
					netCDF::NcVarAtt ua=time_var.getAtt("units");
					if (!ua.isNull())
					{
						units.resize(256,'\0');
						ua.getValues(units.data());
						units.erase(std::find(units.begin(),units.end(),'\0'),units.end());
					}
				}
				catch (...) {}

				auto to_lower=[](std::string s){std::transform(s.begin(),s.end(),s.begin(),::tolower);return s;};
				std::string ul=to_lower(units);
				double mul=1.;
				if (ul.find("hour")!=std::string::npos)     mul=3600.;
				else if (ul.find("day")!=std::string::npos) mul=86400.;
				else if (ul.find("sec")!=std::string::npos) mul=1.;
				// (If other bases, caller must ensure consistency of the epoch.)
				for (auto &t : times) t=static_cast<TimeType>(std::llround(t*mul));  // quantize to integer seconds

				time_values_[file_path]=times; // Store times with file_path as key
				all_time_values_.insert(times.begin(),times.end());
			}
			else
			{
				std::cerr<<"Time variable not found in file: "<<file_path<<std::endl;
			}
		}
		catch (const std::exception& e)
		{
			std::cerr<<"Error reading time variable from file '"<<file_path<<"': "<<e.what()<<std::endl;
		}

		// Iterate over the variable names to load
		for (const auto& var_name : variable_names_)
		{
			try
			{
				netCDF::NcVar var=dataFile.getVar(var_name);
				if (var.isNull())
				{
					throw std::runtime_error("Variable not found: "+var_name);
				}

				// Get variable dimensions
				std::vector<std::size_t> shape;
				[[maybe_unused]] std::size_t size=1;
				for (int i=0;i<var.getDimCount();i++)
				{
					std::size_t dim_size=var.getDim(i).getSize();
					shape.push_back(dim_size);
					size*=dim_size;
				}

				// Read variable attributes for scale_factor and add_offset
				double scale_factor=1.;
				double add_offset=0.;

				// Attempt to read 'scale_factor' attribute
				try
				{
					netCDF::NcVarAtt scale_factor_att=var.getAtt("scale_factor");
					scale_factor_att.getValues(&scale_factor);
				} catch (netCDF::exceptions::NcException& e)
				{
					// Attribute not found; use default scale_factor = 1.0
					scale_factor=1.;
				}

				// Attempt to read 'add_offset' attribute
				try
				{
					netCDF::NcVarAtt add_offset_att=var.getAtt("add_offset");
					add_offset_att.getValues(&add_offset);
				} catch (netCDF::exceptions::NcException& e)
				{
					// Attribute not found; use default add_offset = 0.0
					add_offset=0.;
				}

				// Read data based on variable type
				netCDF::NcType varType=var.getType();
				std::string typeName=varType.getName();

				// Check if the variable type matches DataType or can be converted to DataType
				if (canConvert(varType))
				{
					readVariableData(var,var_name,shape,scale_factor,add_offset,file_path);
				}
				else
				{
					std::cerr<<"Variable '"<<var_name<<"' is of type '"<<typeName<<"', which does not match template type."<<std::endl;
					continue;
				}

				// Store variable shape and type
				shapes_[var_name]=shape;
				types_[var_name]=typeName;

			}
			catch (const std::exception& e)
			{
				std::cerr<<"Error loading variable '"<<var_name<<"' from file '"<<file_path<<"': "<<e.what()<<std::endl;
			}
		}
	}

	void readVariableData(netCDF::NcVar& var,const std::string& var_name,const std::vector<std::size_t>& shape,double scale_factor,double add_offset,const std::string& file_path)
	{
		netCDF::NcType varType=var.getType();

		// Get the NetCDF type ID
		auto typeId=varType.getId();

		// Map NetCDF type to C++ type using a switch statement
		switch (typeId)
		{
			case NC_BYTE:
			{
				using VarDataType=signed char;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			case NC_CHAR:
			{
				using VarDataType=char;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			case NC_SHORT:
			{
				using VarDataType=short;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			case NC_INT:
			{
				using VarDataType=int;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			case NC_INT64:
			{
				using VarDataType=long long;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			case NC_FLOAT:
			{
				using VarDataType=float;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			case NC_DOUBLE:
			{
				using VarDataType=double;
				readAndTransformData<VarDataType>(var,var_name,shape,scale_factor,add_offset,file_path);
				break;
			}
			default:
				throw std::runtime_error("Unsupported NetCDF variable type.");
		}
	}

	template<typename VarDataType>
	void readAndTransformData(netCDF::NcVar& var,const std::string& var_name,const std::vector<std::size_t>& shape,double scale_factor,double add_offset,const std::string& file_path)
	{
		// Read missing value as VarDataType
		VarDataType missing_value=VarDataType{};
		bool has_missing_value=false;

		// Attempt to read '_FillValue' or 'missing_value' attribute
		try
		{
			netCDF::NcVarAtt missing_value_att=var.getAtt("_FillValue");
			if (!missing_value_att.isNull())
			{
				missing_value_att.getValues(&missing_value);
				has_missing_value=true;
			}
			else
			{
				missing_value_att=var.getAtt("missing_value");
				if (!missing_value_att.isNull())
				{
					missing_value_att.getValues(&missing_value);
					has_missing_value=true;
				}
			}
		}
		catch (const netCDF::exceptions::NcException&)
		{
			// Attribute not found; proceed without missing value
		}

		// --- Fallback: treat -30000 as missing if no attribute was found ---
		// (Skip obvious coordinate vars.)
		if (!has_missing_value && var_name!="lon" && var_name!="lat" && var_name!="depth")
		{
			missing_value=static_cast<VarDataType>(-30000);
			has_missing_value=true;
		}

		// Determine if the variable has a time dimension
		int time_dim_index=-1;
		for (int i=0;i<var.getDimCount();i++)
		{
			if (var.getDim(i).getName()=="time")
			{
				time_dim_index=i;
				break;
			}
		}

		if (time_dim_index!=-1)
		{
			// Variable has a time dimension
			// Get the time values for this file
			const auto& times=time_values_[file_path];
			std::size_t time_size=shape[time_dim_index];

			// Prepare start and count vectors for reading data
			std::vector<std::size_t> start(var.getDimCount(),0);
			std::vector<std::size_t> count=shape; // Copy shape
			count[time_dim_index]=1; // We will read one time step at a time

			// Loop over each time index
			for (std::size_t t=0;t<time_size;t++)
			{
				start[time_dim_index]=t;
				TimeType time_value=times[t];

				// Compute the size of data to read for this time step
				std::size_t data_size=1;
				for (std::size_t dim_size : count)
				{
					data_size*=dim_size;
				}

				// Read data for the current time step
				std::vector<VarDataType>raw_data(data_size);
				var.getVar(start,count,raw_data.data());

				// Apply scale_factor and add_offset, handle missing values
				std::vector<DataType> transformed_data(data_size);
				for (std::size_t i=0;i<data_size;i++)
				{
					if (has_missing_value && raw_data[i]==missing_value)
					{
						transformed_data[i]=std::numeric_limits<DataType>::quiet_NaN();
						//transformed_data[i]=std::numeric_limits<DataType>::lowest();
						//std::cerr<<raw_data[i]<<'\t'<<missing_value<<std::endl;
					}
					else
					{
					transformed_data[i]=static_cast<DataType>(raw_data[i])*scale_factor+add_offset;
					}
				}

				// Store the transformed data
				data_[var_name][time_value]=std::move(transformed_data);
			}
		}
		else
		{
			// Variable without a time dimension
			// Compute the total size
			std::size_t size=1;
			for (std::size_t dim_size : shape)
			{
				size*=dim_size;
			}

			// Read data
			std::vector<VarDataType> raw_data(size);
			var.getVar(raw_data.data());

			// Apply scale_factor and add_offset, handle missing values
			std::vector<DataType> transformed_data(size);
			for (std::size_t i=0;i<size;i++)
			{
				if (has_missing_value && raw_data[i]==missing_value)
				{
					transformed_data[i]=std::numeric_limits<DataType>::quiet_NaN();
					//transformed_data[i]=std::numeric_limits<DataType>::lowest();
					//std::cerr<<raw_data[i]<<'\t'<<missing_value<<std::endl;
				}
				else
				{
					transformed_data[i]=static_cast<DataType>(raw_data[i])*scale_factor+add_offset;
				}
			}

			// Normalize longitude if variable is 'lon'
			if (var_name=="lon")
			{
				std::transform(transformed_data.begin(),transformed_data.end(),transformed_data.begin(),&NetCDFDataLoader<DataType,TimeType>::normalizeLongitude);
			}

			// Use a default time_value
			TimeType time_value=TimeType{};

			// Store the transformed data
			data_[var_name][time_value]=std::move(transformed_data);
		}
	}

	bool canConvert(const netCDF::NcType& varType)
	{
		// Allow conversion if DataType is arithmetic and varType is numeric
		return std::is_arithmetic<DataType>::value && (
		varType==netCDF::ncByte ||
		varType==netCDF::ncChar ||
		varType==netCDF::ncShort ||
		varType==netCDF::ncInt ||
		varType==netCDF::ncInt64 ||
		varType==netCDF::ncFloat ||
		varType==netCDF::ncDouble
		);
	}

	static DataType normalizeLongitude(DataType lon)
	{
		lon=std::fmod(lon+DataType(180.),DataType(360.));
		if (lon<DataType(0.)) lon+=DataType(360.);
		return lon-DataType(180.);
	}

	std::vector<std::string> file_paths_;
	std::vector<std::string> variable_names_;
	std::map<std::string, std::map<TimeType, std::vector<DataType>>> data_; // Stores *transformed* data_[variable_name][time_value]
	std::map<std::string,std::vector<std::size_t>> shapes_; // Stores variable shapes
	std::map<std::string,std::string> types_; // Stores variable types
	std::map<std::string,std::vector<TimeType>> time_values_; // Stores time values for each file
	std::set<TimeType> all_time_values_; // Set of all unique time values across files
};

template<typename DataType,typename TimeType=double>
class NetCDFProcessor
{
	public:
	using DataMap=std::map<TimeType,std::vector<DataType>>;
	using VariableProcessingFunc=std::function<void(const std::string&,TimeType,const std::vector<DataType>&)>;

	NetCDFProcessor(const std::string& file_path,const std::vector<std::string>& variable_names)
	: file_path_(file_path),variable_names_(variable_names),data_loader_({file_path},variable_names)
	{
		// Load data during construction
		try
		{
		    data_loader_.loadData();
		}
		catch (const std::exception& e)
		{
			std::cerr<<"Error loading data from file '"<<file_path_<<"': "<<e.what()<<std::endl;
			throw; // Rethrow exception if you want to handle it outside
		}
	}

	void processVariables(const VariableProcessingFunc& processFunc)
	{
		try
		{
			// Iterate over each variable
			for (const auto& var_name : variable_names_)
			{
				// Get data for the variable
				const auto& var_data=data_loader_.getVariableData(var_name);

				// Iterate over time values
				for (const auto& [time_value,data_vector] : var_data)
				{
					// Call the processing function
					processFunc(var_name,time_value,data_vector);
				}
			}
		}
		catch (const std::exception& e)
		{
			std::cerr<<"Error processing variables from file '"<<file_path_<<"': "<<e.what()<<std::endl;
			throw; // Rethrow exception if needed
		}
	}

	// Method to get the underlying netCDF variable (needed for reading attributes)
	netCDF::NcVar getNetCDFVariable(const std::string& variable_name) const
	{
		return data_loader_.getNetCDFVariable(variable_name);
	}

	// Method to get variable shape
	const std::vector<std::size_t>& getVariableShape(const std::string& variable_name) const
	{
		return data_loader_.getVariableShape(variable_name);
	}

	// Method to get variable type
	const std::string& getVariableType(const std::string& variable_name) const
	{
		return data_loader_.getVariableType(variable_name);
	}

	// Method to get time values for this file
	const std::vector<TimeType>& getTimeValues() const
	{
		return data_loader_.getTimeValues(file_path_);
	}

	// Method to get variable data
	const typename NetCDFDataLoader<DataType,TimeType>::DataMap& getVariableData(const std::string& variable_name) const
	{
		return data_loader_.getVariableData(variable_name);
	}

	std::vector<DataType> getVariableDataAtTime(const std::string& var_name,TimeType time_value)
	{
		// Assuming you have a method to get the index of the time_value
		size_t time_index=getTimeIndex(time_value);
		// Read data for the variable at the specified time index
		return data_loader_.getVariableDataAtIndex(var_name,time_index);
	}

	private:
	std::string file_path_;
	std::vector<std::string> variable_names_;
	NetCDFDataLoader<DataType,TimeType> data_loader_; // Stored as a member variable

	// Implement a method to find the index of a given time_value
	//size_t getTimeIndex(TimeType time_value)
	//{
	//	const auto& time_values=data_loader_.getTimeValues(file_path_);
	//	auto it=std::find(time_values.begin(),time_values.end(),time_value);
	//	if (it!=time_values.end())
	//	{
	//		return std::distance(time_values.begin(),it);
	//	}
	//	else
	//	{
	//		throw std::runtime_error("Time value not found in data.");
	//	}
	//}

	size_t getTimeIndex(TimeType time_value)
	{
		const auto& time_values=data_loader_.getTimeValues(file_path_);
		const long long q=std::llround(time_value);
		auto it=std::find_if(time_values.begin(),time_values.end(),[q](TimeType t){return std::llround(t)==q;});
		if (it!=time_values.end()) return static_cast<size_t>(std::distance(time_values.begin(),it));
		throw std::runtime_error("Time value not found in data.");
	}

};

namespace ncdf_utils
{

	//---------------------------------------------------------------------
	// Function: checkFileHeaders
	//
	// Opens a netCDF file in read-only mode and checks its header information
	// without loading full data arrays. It verifies that:
	//   - The "time" variable exists and its first (and only) dimension equals
	//     the expected time count (e.g. 8).
	//   - For each required variable (e.g. "water_u", "water_v", "lat", "lon",
	//     "depth"), the variable exists, has at least one dimension, none of
	//     its dimensions are zero, and if the variable is time-dependent (more
	//     than one dimension), its first dimension matches expectedTimeCount.
	//---------------------------------------------------------------------

	inline std::vector<std::string> checkFileHeaders(const std::string &fileName,const std::vector<std::string> &requiredVariables,int expectedTimeCount)
	{
		std::vector<std::string> errors;
		try
		{
		// Open the netCDF file in read-only mode.
		netCDF::NcFile dataFile(fileName,netCDF::NcFile::read);

		// --- Check the "time" variable ---
		netCDF::NcVar timeVar=dataFile.getVar("time");
		if (timeVar.isNull())
		{
			errors.push_back("Missing variable 'time'");
		}
		else
		{
			if (timeVar.getDimCount()<1)
			{
				errors.push_back("Variable 'time' has no dimensions.");
			}
			else
			{
				std::size_t timeSize=timeVar.getDim(0).getSize();
				if (timeSize!=static_cast<std::size_t>(expectedTimeCount))
				{
					errors.push_back("Variable 'time' dimension size is "+std::to_string(timeSize)+" but expected "+std::to_string(expectedTimeCount));
				}
			}
		}

		// --- Check each other required variable ---
		for (const auto &varName : requiredVariables)
		{
			// Skip "time" because it was already checked.
			if (varName=="time")
			continue;
			netCDF::NcVar var=dataFile.getVar(varName);
			if (var.isNull())
			{
				errors.push_back("Missing variable '"+varName+"'");
			}
			else
			{
				int dimCount=var.getDimCount();
				if (dimCount==0)
				{
					errors.push_back("Variable '"+varName+"' has no dimensions.");
				}
				else
				{
					// Ensure no dimension is zero (i.e. no missing rows or columns).
					for (int i=0;i<dimCount;i++)
					{
						std::size_t dSize=var.getDim(i).getSize();
						if (dSize==0)
						{
							errors.push_back("Variable '"+varName+"' has dimension "+std::to_string(i)+" equal to 0.");
						}
					}
					// For time-dependent variables (those with more than one dimension),
					// assume the first dimension is time and verify its size.
					if (dimCount>1)
					{
						std::size_t timeDim=var.getDim(0).getSize();
						if (timeDim!=static_cast<std::size_t>(expectedTimeCount))
						{
							errors.push_back("Variable '"+varName+"' has time dimension "+std::to_string(timeDim)+" but expected "+std::to_string(expectedTimeCount));
						}
					}
				}
			}
		}
		}
		catch (const std::exception &e)
		{
			errors.push_back("Error processing file '"+fileName+"': "+std::string(e.what()));
		}
		return errors;
	}

} // namespace ncdf_utils

#endif // NCDF_UTILS_HPP
