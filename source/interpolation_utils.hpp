/*
 *  This header provides a small framework for:
 *  -------------------------------------------
 *
 *    1. Time interpolation of scalar fields using cubic B-splines.
 *    2. Spatial interpolation on a sphere (lon/lat) using a KD-tree with
 *       a custom great circle distance metric.
 *    3. Experimental estimation of spatial derivatives of vector/scalar
 *       fields using ``neighbors of neighbors''.
 *
 *  Assumptions and conventions
 *  ---------------------------
 *  - Coordinates are in radians: lon ∈ (−pi, pi], lat ∈ [−pi/2, pi/2].
 *  - GCD and GCD_deriv are assumed to operate in radians and must be
 *    provided by geography_utils.hpp.
 *  - ``shape'' is the Gaussian length-scale parameter used in
 *    exp(−dist^2 / shape). It controls how local the interpolation is and is chosen as the squared grid spacing.
 *  - The radius search parameter rad2 is the squared search radius in the
 *    metric defined by LonLatDistanceAdaptor, it is typically set to (1.1^2)*shape in order to catch at least 4 neighbors.
 *  - The derivative estimation is experimental and not guaranteed to be
 *    formally consistent or optimized. It should be treated as a heuristic
 *    and validated against analytical or high-resolution reference data.
 *  - No internal synchronization is provided. If used from multiple threads,
 *    the caller is responsible for protecting shared data structures.
 */

#include <vector>
#include <array>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "oneapi/tbb/concurrent_vector.h"
#include "nanoflann.hpp"
#include "geography_utils.hpp" 		// provides GCD() and GCD-deriv()
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

namespace bm=boost::math;
namespace nf=nanoflann;

// The type for a 2D point
template<typename DataType>
using point=std::array<DataType,2>;

// Time intervals for spline creation
constexpr double dtt=10800.;		// 3 hours -- exactly like in the HYCOM data, so splines are not really needed here

/*  - TimeSpline<DataType,TimeType>
 *      Wrapper around boost::math::cubic_b_spline.
 *      Stores a 1D time series (values) and constructs a spline over it.
 *      Assumes uniform temporal spacing (time_step), e.g. 3-hour HYCOM data
 *      (see dtt).
 */
template<typename DataType,typename TimeType=double>
class TimeSpline
{
	public:
	std::vector<DataType> values;
	bm::cubic_b_spline<DataType> spline;

	TimeSpline()=default;

	// Function to create the spline
	void createSpline(TimeType start_time,TimeType time_step)
	{
		// Ensure that values are populated
		if (values.empty())
		{
			throw std::runtime_error("No data to create spline.");
		}
		spline=bm::cubic_b_spline<DataType>(values.data(),values.size(),start_time,time_step);
	}
};

/*  - PointCloud<DataType>
 *      Flat container for 2D points on the sphere in (lon, lat) form:
 *          pts = [lon0, lat0, lon1, lat1, ..., lonN-1, latN-1]
 *      Provides the interface expected by nanoflann’s KDTreeSingleIndexAdaptor.
 */
template<typename DataType>
struct PointCloud
{
	// Stores points as [lon0,lat0,lon1,lat1,...,lonN-1,latN-1]
	std::vector<DataType> pts;

	inline std::size_t kdtree_get_point_count() const
	{
		return pts.size()/2;
	}

	// Return the coordinate of point idx (0 for lon, 1 for lat)
	inline DataType kdtree_get_pt(const std::size_t idx,const std::size_t dim) const
	{
		return pts[2*idx+dim];
	}

	// Provide direct access to the contiguous data.
	inline const DataType* data() const {return pts.data();}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const
	{
		return false;
	}
};

/*  - LonLatDistanceAdaptor<DataType, DatasetAdaptor>
 *      Custom metric adaptor for nanoflann that:
 *        * Wraps longitude differences across the dateline (−pi, pi].
 *        * Scales longitudinal distance by cos(average latitude),
 *          approximating distances on the sphere in radians.
 *      Returns squared ``planarized'' distance suitable for KD-tree searches,
 *      but actual great-circle distances are computed separately via GCD().
 */
template <typename DataType,class DatasetAdaptor>
struct LonLatDistanceAdaptor
{
	typedef DataType ElementType;
	typedef DataType DistanceType;

	const DatasetAdaptor &data;

	LonLatDistanceAdaptor(const DatasetAdaptor &data_) : data(data_) {}

	// Accessor needed by nanoflann
	inline const DataType* getPoint(size_t idx) const
	{
		return &data.pts[2*idx];  // accessor
	}

	inline DistanceType evalMetric(const DataType* a,const DataType* b,size_t /*size*/) const
	{
		DataType dlon=a[0]-b[0];
		if (dlon<-pi) dlon+=2*pi;
		else if (dlon>pi) dlon-=2*pi;
		DataType dlat=a[1]-b[1];
		DataType avgLat=(a[1]+b[1])*static_cast<DataType>(.5);
		DataType scaled_dlon=dlon*std::cos(avgLat);
		return scaled_dlon*scaled_dlon+dlat*dlat;
	}

	inline DistanceType evalMetric(const DataType* a,size_t idx_b,size_t size) const
	{
		const DataType* b=getPoint(idx_b);
		return evalMetric(a,b,size);
	}

	template<typename T,typename IndexType>
	inline DistanceType accum_dist(const T a,const T b,IndexType dim) const
	{
		DistanceType diff=static_cast<DistanceType>(a-b);
		if (dim==0)
		{
			// longitude
			if (diff<-pi) diff+=2*pi;
			else if (diff>pi) diff-=2*pi;
		}
		return diff*diff;
	}
};

/*  - KDTree<DataType>
 *      Convenience alias for nanoflann’s KDTreeSingleIndexAdaptor using
 *      the LonLatDistanceAdaptor and PointCloud. This is the spatial index
 *      used for neighbor queries.
 */
template<typename DataType>
using KDTree=nf::KDTreeSingleIndexAdaptor<LonLatDistanceAdaptor<DataType,PointCloud<DataType>>,PointCloud<DataType>,2>;

/*  - Neighbor<DataType>
 *      Container for a single neighbor on the grid:
 *        * ind    : index on the grid (xx_sparse, yy_sparse)
 *        * dist   : great-circle distance (from GCD())
 *        * weight : Gaussian weight exp(−dist^2 / shape)
 */
template<typename DataType>
struct Neighbor
{
	std::size_t ind{};
	DataType dist{};	// great circle distance from GCD function
	DataType weight{};	// weight of the node
	Neighbor(std::size_t index=0,DataType distance=0.,DataType node_weight=0.) : ind(index),dist(distance),weight(node_weight) {}
};

/*  - NeighborsData<DataType, TimeType>
 *      Interface for:
 *        * Querying spatial neighbors within a given search radius (rad2 - squared).
 *        * Computing Gaussian weights based on great-circle distance.
 *        * Interpolating multiple variables in space.
 */
template<typename DataType,typename TimeType=double>
class NeighborsData
{
	private:
	const KDTree<DataType>& index;
	const std::vector<DataType>& xx_sparse{};
	const std::vector<DataType>& yy_sparse{};
	const DataType rad2{},shape{};

	//std::vector<nanoflann::ResultItem<unsigned int,DataType>> indices_dists;

	public:
	NeighborsData(const KDTree<DataType>& idx,const std::vector<DataType>& xxs,const std::vector<DataType>& yys,DataType r2,DataType shape_param)
	: index(idx),xx_sparse(xxs),yy_sparse(yys),rad2(r2),shape(shape_param)
	{}

	 /* Finds neighbors around (x, y) using the KD-tree, computes
	 * great-circle distances and Gaussian weights.
	 */
	void computeNeighbors(std::vector<Neighbor<DataType>>& neighbors,DataType x,DataType y)
	{
		while (x<=-pi) x+=2*pi;
		while (x>  pi) x-=2*pi;
		nanoflann::SearchParameters searchparams;
		searchparams.sorted=false;	// don't sort the results by distance
		//indices_dists.clear();
		std::vector<nanoflann::ResultItem<unsigned int,DataType>> indices_dists;
		neighbors.clear();
		point<DataType> query_point={x,y};
		index.radiusSearch(query_point.data(),rad2,indices_dists,searchparams);
		neighbors.reserve(indices_dists.size());	// Reserve memory to reduce reallocations
		for (const auto& result : indices_dists)
		{
			Neighbor<DataType> nbr;
			nbr.ind=result.first;
			// Check bounds for xx_sparse and yy_sparse
			//if (nbr.ind>=xx_sparse.size() || nbr.ind>=yy_sparse.size())
			//{
			//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in compute method.");
			//}
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=exp(-nbr.dist*nbr.dist/shape);		// Gaussian weighting
			neighbors.push_back(nbr);
		}
	}

	/* Reuses an existing neighbor list but recomputes weights for
         *   a new query position (x, y). Neighbor indices are unchanged.
	 */
	void updateNeighborWeights(std::vector<Neighbor<DataType>>& neighbors,DataType x,DataType y)
	{
		point<DataType> query_point={x,y};
		for (auto& nbr : neighbors)
		{
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=std::exp(-nbr.dist*nbr.dist/shape);
		}
	}

	/* Given:
         * - time t,
         * - precomputed neighbors,
         * - a map var_name -> vector<TimeSpline> (one spline per grid point),
         * it evaluates each spline at time t, computes a weighted average
         * over all neighbors, and normalizes by the total weight.
	 */
	std::map<std::string,DataType> interpolateVariables(
		TimeType t,
		const std::vector<Neighbor<DataType>>& neighbors,
		const std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
		const std::vector<std::string>& variable_names)
	{
		std::map<std::string,DataType> interpolated_values;
		DataType den=static_cast<DataType>(0);

		// Initialize accumulators for each variable
		for (const auto& var_name : variable_names)
		{
			interpolated_values[var_name]=static_cast<DataType>(0);
		}

		for (const auto& nbr : neighbors)
		{
			for (const auto& var_name : variable_names)
			{
				const auto& splines=variable_splines.at(var_name);

				// Check bounds for splines
				if (nbr.ind>=splines.size())
				{
					throw std::out_of_range("Index "+std::to_string(nbr.ind)+
					" out of bounds for variable '"+var_name+
					"' in interpolateVariables method.");
				}

				// Evaluate the spline at time t
				DataType value=splines[nbr.ind].spline(t);

				// Accumulate weighted values
				interpolated_values[var_name]+=nbr.weight*value;
			}

			// Accumulate the total weight
			den+=nbr.weight;
		}

		// Normalize the interpolated values
		if (den>static_cast<DataType>(1e-8))
		{
			for (const auto& var_name : variable_names)
			{
				interpolated_values[var_name]=interpolated_values[var_name]/den;
			}
		}
		else
		{
			// Handle zero denominator case if needed
			for (const auto& var_name : variable_names)
			{
				interpolated_values[var_name]=static_cast<DataType>(0);
			}
		}

		return interpolated_values;
	}

	/* Convenience wrappers combining neighbor computation / update
         *   with interpolation for a given time and position.
	 */
	std::map<std::string,DataType> calculateVariables(
		TimeType t,
		std::vector<Neighbor<DataType>>& neighbors,
		const std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
		const std::vector<std::string>& variable_names,
		DataType x,DataType y)
	{
		computeNeighbors(neighbors,x,y);
		return interpolateVariables(t,neighbors,variable_splines,variable_names);
	}

	std::map<std::string,DataType> calculateVariablesUpdate(
		TimeType t,
		std::vector<Neighbor<DataType>>& neighbors,
		const std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
		const std::vector<std::string>& variable_names,
		DataType x,DataType y)
	{
		updateNeighborWeights(neighbors,x,y);
		return interpolateVariables(t,neighbors,variable_splines,variable_names);
	}
};

// -----------------------------------------------------------------------------
// Experimental attempt to estimate the spatial derivatives of the vector fields
// -----------------------------------------------------------------------------

/*  - NeighborPrime<DataType>
 *      ``Neighbor of a neighbor'':
 *        * ind_prime    : index of the secondary neighbor
 *        * dist_prime   : array [d, d/dx, d/dy] from GCD_deriv()
 *        * weight_prime : Gaussian weight based on d
 */
template<typename DataType>
struct NeighborPrime
{
	std::size_t ind_prime{};
	std::array<DataType,3> dist_prime{};	// great circle distance from GCD function and its derivatives in x and y directions
	DataType weight_prime{};		// weight of the node
	NeighborPrime(std::size_t index_prime=0,std::array<DataType,3> distance_prime={0.,0.,0.},DataType node_weight_prime=0.)
	: ind_prime(index_prime),dist_prime(distance_prime),weight_prime(node_weight_prime) {}
};

/*  - NeighborWithNeighbors<DataType>
 *      Extends Neighbor with a vector<NeighborPrime>:
 *        * Represents one primary neighbor plus its local neighborhood.
 *        * Used to approximate spatial derivatives as a double sum
 *          over neighbors and neighbors-of-neighbors.
 */
template<typename DataType>
struct NeighborWithNeighbors : Neighbor<DataType>		// Derivative is approximated as a double sum
{
	std::vector<NeighborPrime<DataType>> neighbors_prime;
	NeighborWithNeighbors()=default;
	NeighborWithNeighbors(std::size_t index,DataType distance,DataType node_weight,const std::vector<NeighborPrime<DataType>>& nbr_prime)
        : Neighbor<DataType>(index,distance,node_weight),neighbors_prime(nbr_prime) {}
};

/*  - VariableDerivatives<DataType>
 *      Bundles for each variable:
 *        * value                : interpolated value
 *        * time_derivative      : d(value)/dt from spline derivative
 *        * spatial_derivative_x : approximate d(value)/dx
 *        * spatial_derivative_y : approximate d(value)/dy
 */
template<typename DataType>
struct VariableDerivatives
{
	DataType value;                   // Interpolated value
	DataType time_derivative;         // Time derivative
	DataType spatial_derivative_x;    // Spatial derivative in x
	DataType spatial_derivative_y;    // Spatial derivative in y
};

/*  - NeighborWithNeighborsData<DataType, TimeType>
 *      Analogous to NeighborsData, but:
 *        * computeNeighbors(...) builds NeighborWithNeighbors elements:
 *          for each primary neighbor, it performs a second radius search
 *          around the neighbor’s location and constructs NeighborPrime objects.
 */
template<typename DataType,typename TimeType=double>
class NeighborWithNeighborsData
{
	private:
	const KDTree<DataType>& index;
	const std::vector<DataType>& xx_sparse{};
	const std::vector<DataType>& yy_sparse{};
	const DataType rad2{},shape{};

	//std::vector<nanoflann::ResultItem<unsigned int,DataType>> indices_dists;

	public:
	NeighborWithNeighborsData(const KDTree<DataType>& idx,const std::vector<DataType>& xxs,const std::vector<DataType>& yys,DataType r2,DataType shape_param)
	: index(idx),xx_sparse(xxs),yy_sparse(yys),rad2(r2),shape(shape_param) {}

	void computeNeighbors(std::vector<NeighborWithNeighbors<DataType>>& neighbors,DataType x,DataType y)
	{
		point<DataType> query_point={x,y};
		nanoflann::SearchParameters searchparams;
		searchparams.sorted=false;	// don't sort the results by distance
		//indices_dists.clear();
		std::vector<nanoflann::ResultItem<unsigned int,DataType>> indices_dists;
		index.radiusSearch(query_point.data(),rad2,indices_dists,searchparams);
		neighbors.clear();
		neighbors.reserve(indices_dists.size());	// Reserve memory to reduce reallocations
		for (const auto& result : indices_dists)
		{
			NeighborWithNeighbors<DataType> nbr;
			nbr.ind=result.first;
			// Check bounds for xx_sparse and yy_sparse
			//if (nbr.ind>=xx_sparse.size() || nbr.ind>=yy_sparse.size())
			//{
			//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in compute method.");
			//}
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=std::exp(-nbr.dist*nbr.dist/shape);

			point<DataType> query_point_prime={xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind)};
			std::vector<nanoflann::ResultItem<unsigned int,DataType>> indices_dists_prime;
			index.radiusSearch(query_point_prime.data(),rad2,indices_dists_prime,searchparams);
			if(indices_dists_prime.empty())
			{
				//std::cout<<x<<'\t'<<y<<'\n';
				continue;
			}
			nbr.neighbors_prime.clear();
			nbr.neighbors_prime.reserve(indices_dists_prime.size());	// Reserve memory to reduce reallocations
			for (const auto& result_prime : indices_dists_prime)
			{
				NeighborPrime<DataType> nbr_prime;
				nbr_prime.ind_prime=result_prime.first;
				// Check bounds for xx_sparse and yy_sparse
				if (nbr_prime.ind_prime>=xx_sparse.size() || nbr_prime.ind_prime>=yy_sparse.size())
				{
					throw std::out_of_range("Index "+std::to_string(nbr_prime.ind_prime)+" out of bounds in compute method.");
				}
				nbr_prime.dist_prime=GCD_deriv(query_point_prime[0],query_point_prime[1],xx_sparse.at(nbr_prime.ind_prime),yy_sparse.at(nbr_prime.ind_prime));
				nbr_prime.weight_prime=std::exp(-nbr_prime.dist_prime[0]*nbr_prime.dist_prime[0]/shape);
				nbr.neighbors_prime.push_back(nbr_prime);
			}
			neighbors.push_back(nbr);
		}
	}

	void updateNeighborWeights(std::vector<NeighborWithNeighbors<DataType>>& neighbors,DataType x,DataType y)
	{
		point<DataType> query_point={x,y};
		for (auto& nbr : neighbors)
		{
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=std::exp(-nbr.dist*nbr.dist/shape);
			// I don't need to recalculate weights of neighbors' neighbors
		}
	}

	/* interpolateDerivatives(...)
	 *   For each variable:
	 *     - Uses spline and spline.prime to get time behavior.
	 *     - Approximates spatial derivatives via weighted differences
	 *       between primary neighbors and their neighbors, combined
	 *       with distance derivatives from GCD_deriv.
	 *     - Applies a (2 / shape) factor for normalization.
	 */
	std::map<std::string,VariableDerivatives<DataType>> interpolateDerivatives(
		TimeType t,
		const std::vector<NeighborWithNeighbors<DataType>>& neighbors,
		const std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
		const std::vector<std::string>& variable_names,
		DataType shape)  // 'shape' is used for normalization
	{
		// Initialize accumulators for each variable
		std::map<std::string,DataType> value_accumulators;
		std::map<std::string,DataType> time_derivative_accumulators;
		std::map<std::string,DataType> spatial_derivative_x_accumulators;
		std::map<std::string,DataType> spatial_derivative_y_accumulators;
		DataType den=static_cast<DataType>(0);

		// Initialize accumulators for each variable
		for (const auto& var_name : variable_names)
		{
			value_accumulators[var_name]=static_cast<DataType>(0);
			time_derivative_accumulators[var_name]=static_cast<DataType>(0);
			spatial_derivative_x_accumulators[var_name]=static_cast<DataType>(0);
			spatial_derivative_y_accumulators[var_name]=static_cast<DataType>(0);
		}

		// Iterate over neighbors
		for (const auto& nbr : neighbors)
		{
			DataType weight=nbr.weight;
			den+=weight;

			// Accumulate values and time derivatives for each variable
			for (const auto& var_name : variable_names)
			{
				const auto& splines=variable_splines.at(var_name);

				// Check bounds for splines
				if (nbr.ind>=splines.size())
				{
					throw std::out_of_range("Index "+std::to_string(nbr.ind)+
					" out of bounds for variable '"+var_name+
					"' in interpolateDerivatives method.");
				}

				// Evaluate spline and its derivative at time t
				DataType value=splines[nbr.ind].spline(t);
				DataType time_derivative=splines[nbr.ind].spline.prime(t);

				// Accumulate weighted values
				value_accumulators[var_name]+=weight*value;
				time_derivative_accumulators[var_name]+=weight*time_derivative;
			}

			// Initialize accumulators for spatial derivatives for this neighbor
			DataType den_prime=static_cast<DataType>(0);
			std::map<std::string,DataType> spatial_derivative_x_prime;
			std::map<std::string,DataType> spatial_derivative_y_prime;

			for (const auto& var_name : variable_names)
			{
				spatial_derivative_x_prime[var_name]=static_cast<DataType>(0);
				spatial_derivative_y_prime[var_name]=static_cast<DataType>(0);
			}

			// Iterate over neighbor's neighbors
			for (const auto& nbr_prime : nbr.neighbors_prime)
			{
				if (nbr_prime.ind_prime!=nbr.ind)
				{
					DataType d_prime=nbr_prime.dist_prime[0];
					DataType dx_prime=nbr_prime.dist_prime[1];
					DataType dy_prime=nbr_prime.dist_prime[2];

					for (const auto& var_name : variable_names)
					{
						const auto& splines=variable_splines.at(var_name);

						// Check bounds for splines
						if (nbr_prime.ind_prime>=splines.size())
						{
							throw std::out_of_range("Index "+std::to_string(nbr_prime.ind_prime)+
							" out of bounds for variable '"+var_name+
							"' in interpolateDerivatives method.");
						}

						// Compute difference in values
						DataType value_nbr=splines[nbr.ind].spline(t);
						DataType value_nbr_prime=splines[nbr_prime.ind_prime].spline(t);
						DataType dv=value_nbr_prime-value_nbr;

						// Accumulate weighted spatial derivatives
						spatial_derivative_x_prime[var_name]+=nbr_prime.weight_prime*(d_prime*dx_prime*dv);
						spatial_derivative_y_prime[var_name]+=nbr_prime.weight_prime*(d_prime*dy_prime*dv);

						den_prime+=nbr_prime.weight_prime;
					}
				}
			}

			// Accumulate spatial derivatives if den_prime is significant
			if (den_prime>static_cast<DataType>(1e-16))
			{
				for (const auto& var_name : variable_names)
				{
					spatial_derivative_x_accumulators[var_name]+=weight*spatial_derivative_x_prime[var_name]/den_prime;
					spatial_derivative_y_accumulators[var_name]+=weight*spatial_derivative_y_prime[var_name]/den_prime;
				}
			}
		}

		// Normalize accumulators and construct results
		std::map<std::string,VariableDerivatives<DataType>> results;

		if (den>static_cast<DataType>(1e-16))
		{
			for (const auto& var_name : variable_names)
			{
				VariableDerivatives<DataType> var_derivatives;

				var_derivatives.value=value_accumulators[var_name]/den;
				var_derivatives.time_derivative=time_derivative_accumulators[var_name]/den;
				var_derivatives.spatial_derivative_x=(static_cast<DataType>(2.)/shape)*(spatial_derivative_x_accumulators[var_name]/den);
				var_derivatives.spatial_derivative_y=(static_cast<DataType>(2.)/shape)*(spatial_derivative_y_accumulators[var_name]/den);

				results[var_name]=var_derivatives;
			}
		}
		else
		{
			// Handle zero denominator case by initializing to zero
			for (const auto& var_name : variable_names)
			{
				results[var_name]=VariableDerivatives<DataType>{};
			}
		}

		return results;
	}

	std::map<std::string,VariableDerivatives<DataType>> calculateDerivatives(
		TimeType t,
		std::vector<NeighborWithNeighbors<DataType>>& neighbors,
		const std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
		const std::vector<std::string>& variable_names,
		DataType x,DataType y,DataType shape)
	{
		// Recompute neighbors based on (x, y)
		computeNeighbors(neighbors,x,y);

		// Interpolate derivatives for the specified variables
		return interpolateDerivatives(t,neighbors,variable_splines,variable_names,shape);
	}

	std::map<std::string,VariableDerivatives<DataType>> calculateDerivativesUpdate(
		TimeType t,
		std::vector<NeighborWithNeighbors<DataType>>& neighbors,
		const std::map<std::string,std::vector<TimeSpline<DataType,TimeType>>>& variable_splines,
		const std::vector<std::string>& variable_names,
		DataType x,DataType y,DataType shape)
	{
		// Update neighbor weights based on the new position (x, y)
		updateNeighborWeights(neighbors,x,y);

		// Interpolate derivatives for the specified variables
		return interpolateDerivatives(t,neighbors,variable_splines,variable_names,shape);
	}

};
