#pragma once

#include <cmath>

#include "LonLat.h"

namespace optics {

struct Vec3 {
    double coords[3];

    Vec3() : coords{0.0, 0.0, 0.0} {}
    Vec3(double x, double y, double z) : coords{x, y, z} {}
    Vec3(double const* xyz) : coords{xyz[0], xyz[1], xyz[2]} {}

    Vec3(LonLat const& p) {
        double sinLon = std::sin(RAD_PER_DEG * p.lon);
        double cosLon = std::cos(RAD_PER_DEG * p.lon);
        double sinLat = std::sin(RAD_PER_DEG * p.lat);
        double cosLat = std::cos(RAD_PER_DEG * p.lat);
        coords[0] = cosLon * cosLat;
        coords[1] = sinLon * cosLat;
        coords[2] = sinLat;
    }

    Vec3(Vec3 const& v) = default;
    Vec3(Vec3&& v) = default;

    Vec3& operator=(Vec3 const& v) = default;
    Vec3& operator=(Vec3&& v) = default;

    double x() const { return coords[0]; }
    double y() const { return coords[1]; }
    double z() const { return coords[2]; }

    Vec3 operator-() const { return Vec3{-x(), -y(), -z()}; }
    Vec3 operator+(Vec3 const& other) const {
        return Vec3{x() + other.x(), y() + other.y(), z() + other.z()};
    }
    Vec3 operator-(Vec3 const& other) const {
        return Vec3{x() - other.x(), y() - other.y(), z() - other.z()};
    }
    Vec3 operator*(Vec3 const& other) const {
        return Vec3{x() * other.x(), y() * other.y(), z() * other.z()};
    }
    Vec3 operator/(Vec3 const& other) const {
        return Vec3{x() / other.x(), y() / other.y(), z() / other.z()};
    }
    Vec3 operator*(double s) const { return Vec3{x() * s, y() * s, z() * s}; }
    Vec3 operator/(double s) const { return Vec3{x() / s, y() / s, z() / s}; }

    double dot(Vec3 const& other) const {
        return x() * other.x() + y() * other.y() + z() * other.z();
    }

    Vec3 cross(Vec3 const& other) const {
        return Vec3{y() * other.z() - z() * other.y(),
                    z() * other.x() - x() * other.z(),
                    x() * other.y() - y() * other.z()};
    }

    inline LonLat lonLat() const;
};

inline LonLat Vec3::lonLat() const {
    double lon = DEG_PER_RAD * std::atan2(y(), x());
    if (lon < 0.0) {
        lon += 360.0;
    }
    double lat = DEG_PER_RAD * std::asin(z());
    if (lat < -90.0) {
        lat = -90.0;
    } else if (lat > 90.0) {
        lat = 90.0;
    }
    return LonLat{lon, lat};
}

inline Vec3 Normalize(Vec3 const& v) {
    double norm = std::sqrt(v.dot(v));
    return v / norm;
}

inline Vec3 EastOf(LonLat const& p) {
    double sinLon = std::sin(RAD_PER_DEG * p.lon);
    double cosLon = std::cos(RAD_PER_DEG * p.lon);
    double sinLat = std::sin(RAD_PER_DEG * p.lat);
    double cosLat = std::cos(RAD_PER_DEG * p.lat);
    return Vec3{-cosLon * sinLat, -sinLon * sinLat, cosLat};
}

inline Vec3 NorthOf(LonLat const& p) {
    double sinLon = std::sin(RAD_PER_DEG * p.lon);
    double cosLon = std::cos(RAD_PER_DEG * p.lon);
    return Vec3{-sinLon, cosLon, 0.0};
}

inline Vec3 Min(Vec3 const& a, Vec3 const& b) {
    return Vec3{std::min(a.x(), b.x()), std::min(a.y(), b.y()), std::min(a.z(), b.z())};
}

inline Vec3 Max(Vec3 const& a, Vec3 const& b) {
    return Vec3{std::max(a.x(), b.x()), std::max(a.y(), b.y()), std::max(a.z(), b.z())};
}

// Returns the square of the euclidian distance between a and b.
inline double SquaredEuclidianDistance(Vec3 const& a, Vec3 const& b) {
    Vec3 v = a - b;
    return v.dot(v);
}

// Returns the squared euclidian distance between two unit vectors separated
// by the given angle in degrees.
inline double SquaredEuclidianDistance(double angle) {
    double d = std::sin(0.5 * angle * RAD_PER_DEG);
    return 4.0 * d * d;
}

// Computes the minimum squared euclidian distance between two unit vectors a, b
// that have their k-th coordinates fixed to s, t.
inline double MinSquaredEuclidianDistance(double s, double t) {
    return 2.0 * (1.0 - s * t - std::sqrt((1.0 - s * s) * (1.0 - t * t)));
}

inline Vec3 operator*(double s, Vec3 const& v) { return v * s; }

inline Vec3& operator+=(Vec3& a, Vec3 const& b) { return a = (a + b); }
inline Vec3& operator-=(Vec3& a, Vec3 const& b) { return a = (a - b); }
inline Vec3& operator*=(Vec3& a, Vec3 const& b) { return a = (a * b); }
inline Vec3& operator/=(Vec3& a, Vec3 const& b) { return a = (a / b); }
inline Vec3& operator*=(Vec3& a, double s) { return a = (a * s); }
inline Vec3& operator/=(Vec3& a, double s) { return a = (a / s); }

inline bool operator==(Vec3 const& a, Vec3 const& b) {
    return a.x() == b.x() && a.y() == b.y() && a.z() == b.z();
}

inline bool operator!=(Vec3 const& a, Vec3 const& b) { return !(a == b); }

}  // namespace optics