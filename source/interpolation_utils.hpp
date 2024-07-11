#include <vector>
#include <array>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "oneapi/tbb/concurrent_vector.h"
#include "nanoflann.hpp"
#include "geography_utils.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

// Namespace aliasing for clarity
namespace bm=boost::math;
namespace nf=nanoflann;

// Define the type for a 2D point
using point=std::array<double,2>;

// Time intervals for spline creation
constexpr double dtt=10800.;

// Class for time-based cubic spline interpolation
struct TimeSpline
{
	std::vector<double> u_values,v_values;
	bm::cubic_b_spline<double> spline_u,spline_v;

	// Constructor, empty as the struct is a POD
	TimeSpline()=default;

	// Creates cubic B-spline interpolators for the data
	void createSpline()
	{
		spline_u=bm::cubic_b_spline<double>(u_values.data(),u_values.size(),0.,dtt);
		spline_v=bm::cubic_b_spline<double>(v_values.data(),v_values.size(),0.,dtt);
	}
};

// Struct to store and manage a collection of 2D points
struct PointCloud
{
	std::vector<point> points;

	// Returns the total number of points
	inline std::size_t kdtree_get_point_count() const
	{
		return points.size();
	}

	// Function required for compatibility with nanoflann, but not implemented here
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const
	{
		return false;
	}

	// Retrieves the value of the specified coordinate of the point at the given index
	inline double kdtree_get_pt(const std::size_t idx,const std::size_t dim) const
	{
		return points[idx][dim];
	}
};

// Define the type for a KD-tree using nanoflann's KDTreeSingleIndexAdaptor for 2D data
using KDTree=nf::KDTreeSingleIndexAdaptor<nf::L2_Simple_Adaptor<double,PointCloud>,PointCloud,2>;

struct Neighbor
{
	std::size_t ind{};
	double dist{};	// great circle distance from GCD function
	double weight{};	// weight of the node
	Neighbor(std::size_t index=0,double distance=0.,double node_weight=0.) : ind(index), dist(distance), weight(node_weight) {}
};

class NeighborsData
{
	private:
	const KDTree& index;
	const std::vector<double>& xx_sparse{};
	const std::vector<double>& yy_sparse{};
	const double rad2{},shape{};

	//std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists;

	public:
	NeighborsData(const KDTree& idx,const std::vector<double>& xxs,const std::vector<double>& yys,double r2,double shape_param)
	: index(idx),xx_sparse(xxs),yy_sparse(yys),rad2(r2),shape(shape_param)
	{}

	void computeNeighbors(std::vector<Neighbor>& neighbors,double x,double y)
	{
		point query_point={x,y};
		nanoflann::SearchParameters searchparams;
		searchparams.sorted=false;	// don't sort the results by distance
		//indices_dists.clear();
		std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists;
		index.radiusSearch(query_point.data(),rad2,indices_dists,searchparams);
		neighbors.clear();
		neighbors.reserve(indices_dists.size());	// Reserve memory to reduce reallocations
		for (const auto& result : indices_dists)
		{
			Neighbor nbr;
			nbr.ind=result.first;
			// Check bounds for xx_sparse and yy_sparse
			//if (nbr.ind>=xx_sparse.size() || nbr.ind>=yy_sparse.size())
			//{
			//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in compute method.");
			//}
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=exp(-nbr.dist*nbr.dist/shape);
			neighbors.push_back(nbr);
		}
		if (x>359./360.*2*M_PI)
		{
			point query_point_margin={x-2*M_PI,y};
			std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists_margin;
			index.radiusSearch(query_point_margin.data(),rad2,indices_dists_margin,searchparams);
			for (const auto& result : indices_dists_margin)
			{
				Neighbor nbr;
				nbr.ind=result.first;
				// Check bounds for xx_sparse and yy_sparse
				//if (nbr.ind>=xx_sparse.size() || nbr.ind>=yy_sparse.size())
				//{
				//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in compute method.");
				//}
				nbr.dist=GCD(query_point_margin[0],query_point_margin[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
				nbr.weight=exp(-nbr.dist*nbr.dist/shape);
				neighbors.push_back(nbr);
			}
		}
		if (x<1./360.*2*M_PI)
		{
			point query_point_margin={x+2*M_PI,y};
			std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists_margin;
			index.radiusSearch(query_point_margin.data(),rad2,indices_dists_margin,searchparams);
			for (const auto& result : indices_dists_margin)
			{
				Neighbor nbr;
				nbr.ind=result.first;
				// Check bounds for xx_sparse and yy_sparse
				//if (nbr.ind>=xx_sparse.size() || nbr.ind>=yy_sparse.size())
				//{
				//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in compute method.");
				//}
				nbr.dist=GCD(query_point_margin[0],query_point_margin[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
				nbr.weight=exp(-nbr.dist*nbr.dist/shape);
				neighbors.push_back(nbr);
			}
		}
	}

	void updateNeighborWeights(std::vector<Neighbor>& neighbors,double x,double y)
	{
		point query_point={x,y};
		for (auto& nbr : neighbors)
		{
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=std::exp(-nbr.dist*nbr.dist/shape);
		}
	}

	std::vector<double> interpolateVelocities(double t,const std::vector<Neighbor>& neighbors,const oneapi::tbb::concurrent_vector<TimeSpline>& splines)
	{
		double uu{0.},vv{0.},den{0.};

		for (const auto& nbr : neighbors)
		{
			// Check bounds for splines
			if (nbr.ind>=splines.size())
			{
				throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in interpolateVelocities method.");
			}
			uu+=nbr.weight*splines.at(nbr.ind).spline_u(t);
			vv+=nbr.weight*splines.at(nbr.ind).spline_v(t);
			den+=nbr.weight;
		}

		std::vector<double> vel(2,0.);
		if (den>1e-8)
		{
			vel.at(0)=uu/den;
			vel.at(1)=vv/den;
		}
		return vel;
	}

	std::vector<double> calculateVelocities(double t,std::vector<Neighbor>& neighbors,const oneapi::tbb::concurrent_vector<TimeSpline>& splines,double x,double y)
	{
		computeNeighbors(neighbors,x,y);
		return interpolateVelocities(t,neighbors,splines);
	}

	std::vector<double> calculateVelocitiesUpdate(double t,std::vector<Neighbor>& neighbors,const oneapi::tbb::concurrent_vector<TimeSpline>& splines,double x,double y)
	{
		updateNeighborWeights(neighbors,x,y);
		return interpolateVelocities(t,neighbors,splines);
	}
};

struct NeighborPrime
{
	std::size_t ind_prime{};
	std::array<double,3> dist_prime{};	// great circle distance from GCD function
	double weight_prime{};	// weight of the node
	NeighborPrime(std::size_t index_prime=0,std::array<double,3> distance_prime={0.,0.,0.},double node_weight_prime=0.)
	: ind_prime(index_prime),dist_prime(distance_prime),weight_prime(node_weight_prime) {}
};

struct NeighborWithNeighbors : Neighbor
{
	std::vector<NeighborPrime> neighbors_prime;
	NeighborWithNeighbors()=default;
	NeighborWithNeighbors(std::size_t index,double distance,double node_weight,const std::vector<NeighborPrime>& nbr_prime)
        : Neighbor(index,distance,node_weight),neighbors_prime(nbr_prime) {}
};

class NeighborWithNeighborsData
{
	private:
	const KDTree& index;
	const std::vector<double>& xx_sparse{};
	const std::vector<double>& yy_sparse{};
	const double rad2{},shape{};

	//std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists;

	public:
	NeighborWithNeighborsData(const KDTree& idx,const std::vector<double>& xxs,const std::vector<double>& yys,double r2,double shape_param)
	: index(idx),xx_sparse(xxs),yy_sparse(yys),rad2(r2),shape(shape_param) {}

	void computeNeighbors(std::vector<NeighborWithNeighbors>& neighbors,double x,double y)
	{
		point query_point={x,y};
		nanoflann::SearchParameters searchparams;
		searchparams.sorted=false;	// don't sort the results by distance
		//indices_dists.clear();
		std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists;
		index.radiusSearch(query_point.data(),rad2,indices_dists,searchparams);
		neighbors.clear();
		neighbors.reserve(indices_dists.size());	// Reserve memory to reduce reallocations
		for (const auto& result : indices_dists)
		{
			NeighborWithNeighbors nbr;
			nbr.ind=result.first;
			// Check bounds for xx_sparse and yy_sparse
			//if (nbr.ind>=xx_sparse.size() || nbr.ind>=yy_sparse.size())
			//{
			//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in compute method.");
			//}
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=std::exp(-nbr.dist*nbr.dist/shape);

			point query_point_prime={xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind)};
			std::vector<nanoflann::ResultItem<unsigned int,double>> indices_dists_prime;
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
				NeighborPrime nbr_prime;
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

	void updateNeighborWeights(std::vector<NeighborWithNeighbors>& neighbors,double x,double y)
	{
	point query_point={x,y};
	for (auto& nbr : neighbors)
		{
			nbr.dist=GCD(query_point[0],query_point[1],xx_sparse.at(nbr.ind),yy_sparse.at(nbr.ind));
			nbr.weight=std::exp(-nbr.dist*nbr.dist/shape);
			// I don't need to recalculate weights of neighbors' neighbors
		}
	}

	std::vector<double> interpolateDerivatives(double t,const std::vector<NeighborWithNeighbors>& neighbors,const oneapi::tbb::concurrent_vector<TimeSpline>& splines)
	{
		double den{0.},uu{0.},vv{0.},dudx{0.},dudy{0.},dvdx{0.},dvdy{0.},dudt{0.},dvdt{0.};
		for (const auto& nbr : neighbors)
		{
			// Check bounds for splines
			//if (nbr.ind>=splines.size())
			//{
			//	throw std::out_of_range("Index "+std::to_string(nbr.ind)+" out of bounds in interpolateDerivatives method.");
			//}
			uu+=nbr.weight*splines.at(nbr.ind).spline_u(t);
			vv+=nbr.weight*splines.at(nbr.ind).spline_v(t);
			uu+=nbr.weight*splines.at(nbr.ind).spline_u.prime(t);
			vv+=nbr.weight*splines.at(nbr.ind).spline_v.prime(t);
			den+=nbr.weight;
			double den_prime{0.},dudx_prime{0.},dudy_prime{0.},dvdx_prime{0.},dvdy_prime{0.};
			for (const auto& nbr_prime : nbr.neighbors_prime)
			{
				//if (std::abs(nbr_prime.dist_prime[0])>1e-16)
				if (nbr_prime.ind_prime!=nbr.ind)
				{
					// Check bounds for splines
					if (nbr_prime.ind_prime>=splines.size())
					{
						throw std::out_of_range("Index "+std::to_string(nbr_prime.ind_prime)+" out of bounds in interpolateDerivatives method.");
					}
					den_prime+=nbr_prime.weight_prime;
					double d_prime=nbr_prime.dist_prime[0];
					double dx_prime=nbr_prime.dist_prime[1],dy_prime=nbr_prime.dist_prime[2];
					double du=splines.at(nbr_prime.ind_prime).spline_u(t)-splines.at(nbr.ind).spline_u(t);
					double dv=splines.at(nbr_prime.ind_prime).spline_v(t)-splines.at(nbr.ind).spline_v(t);
					dudx_prime+=nbr_prime.weight_prime*(d_prime*dx_prime*du);
					dudy_prime+=nbr_prime.weight_prime*(d_prime*dy_prime*du);
					dvdx_prime+=nbr_prime.weight_prime*(d_prime*dx_prime*dv);
					dvdy_prime+=nbr_prime.weight_prime*(d_prime*dy_prime*dv);
					//std::cout<<nbr_prime.weight_prime<<'\t'<<nbr_prime.dist_prime[0]<<'\t'<<nbr_prime.dist_prime[1]<<'\t'<<nbr_prime.dist_prime[2]<<'\t'<<nbr.ind<<'\n';
				}
				else
				{
					dudx_prime+=0.;
					dudy_prime+=0.;
					dvdx_prime+=0.;
					dvdy_prime+=0.;
				}

			}
			if (den_prime>1e-16)
			{
				dudx+=nbr.weight*dudx_prime/den_prime;
				dudy+=nbr.weight*dudy_prime/den_prime;
				dvdx+=nbr.weight*dvdx_prime/den_prime;
				dvdy+=nbr.weight*dvdy_prime/den_prime;
			}
		}

		std::vector<double> results(8,0.);
		if (den>1e-16)
		{
			results[0]=uu/den;
			results[1]=vv/den;
			results[2]=2./shape*dudx/den;
			results[3]=2./shape*dudy/den;
			results[4]=2./shape*dvdx/den;
			results[5]=2./shape*dvdy/den;
			results[6]=dudt/den;
			results[7]=dvdt/den;
		}
		return results;
	}

	std::vector<double> calculateDerivatives(double t,std::vector<NeighborWithNeighbors>& neighbors,const oneapi::tbb::concurrent_vector<TimeSpline>& splines,double x,double y)
	{
		computeNeighbors(neighbors,x,y);
		return interpolateDerivatives(t,neighbors,splines);
	}

	std::vector<double> calculateDerivativesUpdate(double t,std::vector<NeighborWithNeighbors>& neighbors,const oneapi::tbb::concurrent_vector<TimeSpline>& splines,double x,double y)
	{
		updateNeighborWeights(neighbors,x,y);
		return interpolateDerivatives(t,neighbors,splines);
	}
};
