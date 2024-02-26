#pragma once

#include <random>
#include <string_view>

namespace optics {

constexpr double PI = 3.1415926535897932384626433832795;
constexpr double ONE_OVER_PI = 0.318309886183790671537767526745;
constexpr double RAD_PER_DEG = 0.0174532925199432957692369076849;
constexpr double DEG_PER_RAD = 57.2957795130823208767981548141;

// A point on the unit sphere (sky), specified in spherical polar coordinates.
//
// All angles (stored or passed to member functions) are in units of degrees.
struct LonLat {
    double lon;
    double lat;

    static LonLat fromDegrees(double lon, double lat) { return LonLat{lon, lat}; }

    static LonLat fromRadians(double lon, double lat) {
        return LonLat{DEG_PER_RAD * lon, DEG_PER_RAD * lat};
    }

    // Picks a point uniformly at random on the unit sphere.
    static LonLat random(std::mt19937_64& rng);

    // Picks a point uniformly at random in the specified latitude range.
    static LonLat random(std::mt19937_64& rng, double latMin, double latMax);

    // Picks a point uniformly at random in the specified longitude/latitude
    // range.
    static LonLat random(std::mt19937_64& rng, double lonMin, double lonMax,
                         double latMin, double latMax);

    // Creates a point from the longitude / RA and latitude / Dec in the first two
    // values of the given CSV string. Values must be in units of degrees, and are
    // assumed not to be escaped or quoted.
    static LonLat fromCsv(std::string_view csv, char delim);

    // Returns a copy of this point randomly perturbed according to a normal
    // distribution centered on the original point and with a standard deviation of
    // `sigma` degrees.
    LonLat perturb(std::mt19937_64& rng, double sigma) const;

    // Returns a copy of this point randomly perturbed in the direction given by the
    // specified position angle, such that the distance to the original point is
    // normally distributed with a standard deviation of `sigma` degrees.
    LonLat perturb(std::mt19937_64& rng, double sigma, double positionAngle) const;

    // Returns the angle between this point and p, in degrees.
    double distance(LonLat const& p) const;
};

// Returns the width in longitude α of minimal lon/lat bounding boxes for small circles
// with the given radius (in deg) and center latitude (also in deg).
double LongitudeExtent(double radius, double lat);

}  // namespace optics
