#pragma once

#include "vec3.h"

class ray {
  private:
    point3 orig{};
    vec3 dir{};
    double tm{};

  public:
    ray() = default;
    ray(const point3 &origin, const vec3 &direction, double time = 0)
        : orig{origin}, dir{direction}, tm(time) {}

    const point3 &origin() const { return orig; }
    const point3 &direction() const { return dir; }
    const double &time() const { return tm; }
    point3 at(double t) const { return orig + t * dir; }
};
