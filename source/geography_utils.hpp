#include <array>
#include <cmath>

double GCD(const double& lon1,const double& lat1,const double& lon2,const double& lat2)
{
	double cn=std::sin(lat1)*std::sin(lat2)+std::cos(lat1)*std::cos(lat2)*std::cos(lon1-lon2);
	double sn=std::sqrt(1.-cn*cn);
	return std::atan2(sn,cn);
}

// Computes the Great Circle Distance and its partial derivatives with respect to lat1 and lon1.
// Params: lat1, lon1 - latitude and longitude of the first point in radians.
//         lat2, lon2 - latitude and longitude of the second point in radians.
// Returns: std::array containing [central angle, partial derivative w.r.t. lat1, partial derivative w.r.t. lon1]
std::array<double,3> GCD_deriv(const double& lon1,const double& lat1,const double& lon2,const double& lat2)
{
	double cn=std::sin(lat1)*std::sin(lat2)+std::cos(lat1)*std::cos(lat2)*std::cos(lon1-lon2);
	double sn=std::sqrt(1.-cn*cn);
	double tn=std::atan2(sn,cn);
	double dx=std::cos(lat1)*std::cos(lat2)*std::sin(lon1-lon2)/sn;
	double dy=(std::sin(lat1)*std::cos(lat2)*std::cos(lon1-lon2)-std::cos(lat1)*std::sin(lat2))/sn;
	return {tn,dx,dy};
}

// Given an offset (x,y) from a reference point with latitude and longitude, this function computes the new geographic coordinates after applying the offset.
// Params: x, y - offsets in meters.
//         lat, lon - reference latitude and longitude in radians. These are modified in-place.
void inverseTransform(const double& x,const double& y,double& lon,double& lat)
{
    constexpr double EARTH_RADIUS=6378000.; // Radius of the Earth in meters.

    double nrm=std::sqrt(x*x+y*y);
    double gcd=nrm/EARTH_RADIUS;
    double cn=std::cos(gcd);
    double sn=std::sin(gcd);
    double cl=std::cos(lat);
    double sl=std::sin(lat);
    double lat_prime=std::asin(sl*cn+y/nrm*cl*sn);
    double lon_prime=lon+std::atan2(x*sn,nrm*cl*cn-y*sl*sn);

    lat=lat_prime;
    lon=lon_prime;
}
