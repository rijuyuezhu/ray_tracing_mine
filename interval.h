#pragma once

#include "utils.h"

#include <algorithm>

class interval {
  public:
    double min;
    double max;

  public:
    interval(double min, double max) : min{min}, max{max} {}
    interval() : min{+infinity}, max{-infinity} {}
    interval(const interval &a, const interval &b) {
        min = std::min(a.min, b.min);
        max = std::max(a.max, b.max);
    }

    double size() const { return max - min; }

    bool contains(double x) const { return min <= x && x <= max; }

    bool surrounds(double x) const { return min <= x && x <= max; }

    double clamp(double x) const { return std::clamp(x, min, max); }

    interval expand(double delta) const {
        auto padding{delta / 2};
        return {min - padding, max + padding};
    }

    static const interval empty;
    static const interval universe;
};

inline const interval interval::empty{+infinity, -infinity};
inline const interval interval::universe{-infinity, +infinity};

inline interval operator+(const interval &ival, double displacement) {
    return interval(ival.min + displacement, ival.max + displacement);
}

inline interval operator+(double displacement, const interval &ival) {
    return ival + displacement;
}
