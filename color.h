#pragma once

#include <iostream>

#include "interval.h"
#include "vec3.h"

using color = vec3;

inline double linear_to_gamma(double linear_component) {
    if (linear_component > 0) {
        return std::sqrt(linear_component);
    } else {
        return 0;
    }
}

inline void write_color(std::ostream &out, const color &pixel_color) {
    auto r{pixel_color.x()};
    auto g{pixel_color.y()};
    auto b{pixel_color.z()};

    r = linear_to_gamma(r);
    g = linear_to_gamma(g);
    b = linear_to_gamma(b);

    static const interval intensity(0.000, 0.999);
    auto rbyte{static_cast<int>(256 * intensity.clamp(r))};
    auto gbyte{static_cast<int>(256 * intensity.clamp(g))};
    auto bbyte{static_cast<int>(256 * intensity.clamp(b))};

    out << rbyte << ' ' << gbyte << ' ' << bbyte << '\n';
}
