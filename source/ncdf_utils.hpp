#include <iostream>
#include <vector>
#include <string>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

template<typename T>
class NetCDFWrapper
{
	private:
	NcFile dataFile;

	NcVar getVar(const std::string& varName)
	{
		NcVar var=dataFile.getVar(varName);
		if (var.isNull())
		{
			throw std::runtime_error("Variable not found: "+varName);
		}
		return var;
	}

	public:
	NetCDFWrapper(const std::string& filename):dataFile(filename,NcFile::read){}

	std::vector<T> getVariableData(const std::string& varName)
	{
		NcVar var=getVar(varName);

		std::size_t size{1};
		for (int i=0;i<var.getDimCount();i++)
		{
			size*=var.getDim(i).getSize();
		}

		std::vector<T> data(size);
		var.getVar(data.data());

		return data;
	}

	std::vector<std::size_t> getVariableShape(const std::string& varName)
	{
		NcVar var=getVar(varName);

		std::vector<std::size_t> shape;
		for (int i=0;i<var.getDimCount();i++)
		{
			shape.push_back(var.getDim(i).getSize());
		}

		return shape;
	}

	std::string getVariableType(const std::string& varName)
	{
		NcVar var=getVar(varName);
		return var.getType().getName();
	}

};
