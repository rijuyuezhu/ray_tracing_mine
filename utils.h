#pragma once

#include <limits>
#include <numbers>
#include <random>

inline constexpr double infinity{std::numeric_limits<double>::infinity()};
using std::numbers::pi;

inline constexpr double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

inline double random_double() {
    thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    thread_local std::mt19937 generator{std::random_device{}()};
    return distribution(generator);
}

inline double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

inline int random_int(int min, int max) {
    return static_cast<int>(random_double(min, max + 1));
}
