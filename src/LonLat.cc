#include "LonLat.h"

#include <fast_float/fast_float.h>
#include <fmt/core.h>

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string_view>

#include "Vec3.h"

namespace optics {

namespace {

void checkLat(double lat) {
    if (lat < -90.0 || lat > 90.0) {
        throw std::invalid_argument(fmt::format("invalid latitude {}", lat));
    }
}

void checkLon(double lon) {
    if (lon < 0.0 || lon > 360.0) {
        throw std::invalid_argument(fmt::format("invalid longitude {}", lon));
    }
}

double randomLat(std::mt19937_64& rng, double latMin, double latMax) {
    checkLat(latMin);
    checkLat(latMax);
    if (latMin > latMax) {
        throw std::invalid_argument(
            fmt::format("invalid latitude range [{}, {}]", latMin, latMax));
    }
    double zMin = std::sin(RAD_PER_DEG * latMin);
    double zMax = std::sin(RAD_PER_DEG * latMax);
    double z = zMin == zMax ? zMin : std::uniform_real_distribution{zMin, zMax}(rng);
    double lat = DEG_PER_RAD * std::asin(z);
    return lat < latMin ? latMin : lat > latMax ? latMax : lat;
}

}  // namespace

LonLat LonLat::random(std::mt19937_64& rng) {
    double lon = std::uniform_real_distribution{0.0, 360.0}(rng);
    double z = std::uniform_real_distribution{-1.0, 1.0}(rng);
    return LonLat{lon, DEG_PER_RAD * std::asin(z)};
}

LonLat LonLat::random(std::mt19937_64& rng, double latMin, double latMax) {
    double lon = std::uniform_real_distribution{0.0, 360.0}(rng);
    double lat = randomLat(rng, latMin, latMax);
    return LonLat{lon, lat};
}

LonLat LonLat::random(std::mt19937_64& rng, double lonMin, double lonMax, double latMin,
                      double latMax) {
    checkLon(lonMin);
    checkLon(lonMax);
    double lat = randomLat(rng, latMin, latMax);
    double lon;
    if (lonMin < lonMax) {
        lon = std::uniform_real_distribution{lonMin, lonMax}(rng);
    } else if (lonMin == lonMax) {
        lon = lonMin;
    } else {
        lon = std::uniform_real_distribution{lonMin - 360.0, lonMax}(rng);
        if (lon < 0.0) {
            lon += 360.0;
        }
    }
    return LonLat{lon, lat};
}

LonLat LonLat::fromCsv(std::string_view csv, char delim) {
    double lon;
    double lat;
    char delimOrNewline[2] = {delim, '\n'};

    auto firstDelim = csv.find_first_of(delim, 0);
    if (firstDelim == std::string_view::npos) {
        throw std::invalid_argument(fmt::format(
            "csv line {} (delim={}) does not begin with lon,lat fields", csv, delim));
    }
    auto lonEnd = csv.begin() + firstDelim;
    auto lonResult = fast_float::from_chars(csv.begin(), lonEnd, lon);
    if (lonResult.ec != std::errc{} || lonResult.ptr != lonEnd || lon < -360.0 ||
        lon > 360.0) {
        throw std::invalid_argument(fmt::format(
            "first field of csv line {} (delim={}) is not a valid longitude", csv,
            delim));
    }
    if (lon < 0.0) {
        lon += 360.0;
    }

    auto secondDelimOrNewline = csv.find_first_of(
        std::string_view{delimOrNewline, sizeof(delimOrNewline)}, firstDelim + 1);
    auto latEnd = (secondDelimOrNewline == std::string_view::npos)
                      ? csv.end()
                      : csv.begin() + secondDelimOrNewline;
    auto latResult = fast_float::from_chars(lonEnd + 1, latEnd, lat);
    if (latResult.ec != std::errc{} || latResult.ptr != lonEnd || lat < -90.0 ||
        lat > 90.0) {
        throw std::invalid_argument(fmt::format(
            "second field of csv line {} (delim={}) is not a valid latitude", csv,
            delim));
    }

    return LonLat{lon, lat};
}

LonLat LonLat::perturb(std::mt19937_64& rng, double sigma) const {
    return perturb(rng, sigma, std::uniform_real_distribution{0.0, 360.0}(rng));
}

LonLat LonLat::perturb(std::mt19937_64& rng, double sigma, double positionAngle) const {
    Vec3 v = *this;
    Vec3 n = NorthOf(*this);
    Vec3 e = EastOf(*this);

    // rotate north vector n at v by -pa
    Vec3 t = std::sin(RAD_PER_DEG * positionAngle) * e +
             std::cos(RAD_PER_DEG * positionAngle) * n;

    // perturb in this direction by a random angle that is normally distributed with a
    // standard deviation of sigma degrees
    double mag = std::normal_distribution{0.0, RAD_PER_DEG * sigma}(rng);
    Vec3 p = std::cos(mag) * v + std::sin(mag) * t;

    return p.lonLat();
}

double LonLat::distance(LonLat const& p) const {
    Vec3 v0(*this);
    Vec3 v1(p);
    return DEG_PER_RAD * std::acos(v0.dot(v1));
}

double LongitudeExtent(double radius, double lat) {
    constexpr double POLE_EPSILON = 1e-6;

    if (radius < 0.0 || radius > 90.0) {
        throw std::invalid_argument("radius must be in range [0, 90.0] deg");
    }
    if (radius == 0.0) {
        return 0.0;
    }
    lat = std::max(-90.0, lat);
    lat = std::min(90.0, lat);
    if (std::fabs(lat) + radius > 90.0 - POLE_EPSILON) {
        return 360.0;
    }
    double y = std::sin(RAD_PER_DEG * radius);
    double x = std::sqrt(std::fabs(std::cos(RAD_PER_DEG * (lat - radius)) *
                                   std::cos(RAD_PER_DEG * (lat + radius))));
    return 2.0 * DEG_PER_RAD * std::fabs(std::atan(y / x));
}

}  // namespace optics
