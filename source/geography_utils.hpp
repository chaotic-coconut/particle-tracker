#include <array>
#include <cmath>

// Utility: clamp a value between lo and hi.
template<typename DataType>
DataType clamp(DataType val,DataType lo,DataType hi)
{
	return std::max(lo,std::min(val,hi));
}

// Computes the great circle distance (central angle) between two points on a sphere,
// given their longitudes and latitudes (in radians).
//
// The formula uses the spherical law of cosines to compute the cosine of the central angle,
// then derives the sine of the angle and returns the angle via atan2 for better numerical stability.
template<typename DataType>
DataType GCD(const DataType& lon1,const DataType& lat1,const DataType& lon2,const DataType& lat2)
{
	DataType cn=std::sin(lat1)*std::sin(lat2)+std::cos(lat1)*std::cos(lat2)*std::cos(lon1-lon2);
	cn=clamp(cn,static_cast<DataType>(-1),static_cast<DataType>(1));
	DataType sn=std::sqrt(1.-cn*cn);
	return std::atan2(sn,cn);
}

// Computes the great circle distance (central angle) and its partial derivatives with respect to the first point's
// latitude (lat1) and longitude (lon1). All angles are in radians.
// Returns an array containing three values:
//   [0] - The central angle (great circle distance) between the two points,
//   [1] - The partial derivative of the central angle with respect to lon1,
//   [2] - The partial derivative of the central angle with respect to lat1.
template<typename DataType>
std::array<DataType,3> GCD_deriv(const DataType& lon1,const DataType& lat1,const DataType& lon2,const DataType& lat2)
{
	DataType cn=std::sin(lat1)*std::sin(lat2)+std::cos(lat1)*std::cos(lat2)*std::cos(lon1-lon2);
	cn=clamp(cn,static_cast<DataType>(-1),static_cast<DataType>(1));
	DataType sn=std::sqrt(1.-cn*cn);
	DataType tn=std::atan2(sn,cn);
	// Compute the partial derivative with respect to longitude (lon1).
	// This derivative quantifies how a small change in lon1 affects the central angle.
	DataType dx=tn*std::cos(lat2)*std::sin(lon2-lon1)/sn;
	// Compute the partial derivative with respect to latitude (lat1).
	// This derivative quantifies how a small change in lat1 affects the central angle.
	DataType dy=tn*(std::cos(lat1)*std::sin(lat2)-std::sin(lat1)*std::cos(lat2)*std::cos(lon2-lon1))/sn;
	return {tn,dx,dy};
}

// Given an offset (x,y) from a reference point with latitude and longitude, this function computes the new geographic coordinates after applying the offset.
// Params: x, y - offsets in meters.
//         lat, lon - reference latitude and longitude in radians. These are modified in-place.
template<typename DataType>
void inverseTransform(const DataType& x,const DataType& y,DataType& lon,DataType& lat)
{
	constexpr DataType EARTH_RADIUS=6378000.; // Radius of the Earth in meters.

	DataType nrm=std::sqrt(x*x+y*y);
	const DataType eps=1.e-10;
	// If the displacement is negligible, do nothing.
	if (nrm<eps) return;

	DataType gcd=nrm/EARTH_RADIUS;
	DataType cn=std::cos(gcd);
	DataType sn=std::sin(gcd);
	DataType cl=std::cos(lat);
	DataType sl=std::sin(lat);
	DataType arg=sl*cn+(y/nrm)*cl*sn;
	arg=clamp(arg,static_cast<DataType>(-1),static_cast<DataType>(1));
	DataType lat_prime=std::asin(arg);
	DataType lon_prime=lon+std::atan2(x*sn,nrm*cl*cn-y*sl*sn);

	lat=lat_prime;
	lon=lon_prime;
}
