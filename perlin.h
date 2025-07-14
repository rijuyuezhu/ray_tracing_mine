#pragma once

#include "utils.h"
#include "vec3.h"
class perlin {
  private:
    static constexpr int POINT_COUNT{256};
    vec3 randvec[POINT_COUNT];
    int perm_x[POINT_COUNT];
    int perm_y[POINT_COUNT];
    int perm_z[POINT_COUNT];

  public:
    perlin() {
        for (int i{0}; i < POINT_COUNT; i++) {
            randvec[i] = unit_vector(vec3::random(-1, 1));
        }
        perlin_generate_perm(perm_x);
        perlin_generate_perm(perm_y);
        perlin_generate_perm(perm_z);
    }

    double noise(const point3 &p) const {
        auto u{p.x() - std::floor(p.x())};
        auto v{p.y() - std::floor(p.y())};
        auto w{p.z() - std::floor(p.z())};

        auto i{static_cast<int>(std::floor(p.x()))};
        auto j{static_cast<int>(std::floor(p.y()))};
        auto k{static_cast<int>(std::floor(p.z()))};

        vec3 c[2][2][2];

        for (int di{0}; di < 2; di++)
            for (int dj{0}; dj < 2; dj++)
                for (int dk{0}; dk < 2; dk++) {
                    c[di][dj][dk] = randvec[perm_x[(i + di) & 255] ^
                                            perm_y[(j + dj) & 255] ^
                                            perm_z[(k + dk) & 255]];
                }
        return perlin_interp(c, u, v, w);
    }
    double turb(const point3 &p, int depth) const {
        auto accum = 0.0;
        auto temp_p = p;
        auto weight = 1.0;

        for (int i = 0; i < depth; i++) {
            accum += weight * noise(temp_p);
            weight *= 0.5;
            temp_p *= 2;
        }

        return std::fabs(accum);
    }

  private:
    static void perlin_generate_perm(int *p) {
        for (int i{0}; i < POINT_COUNT; i++) {
            p[i] = i;
        }
        permute(p, POINT_COUNT);
    }

    static void permute(int *p, int n) {
        for (int i{n - 1}; i > 0; i--) {
            int target{random_int(0, i)};
            std::swap(p[i], p[target]);
        }
    }

    static double perlin_interp(const vec3 c[2][2][2], double u, double v,
                                double w) {
        double uu{u * u * (3 - 2 * u)};
        double vv{v * v * (3 - 2 * v)};
        double ww{w * w * (3 - 2 * w)};
        double accum{0.0};

        for (int i{0}; i < 2; i++) {
            for (int j{0}; j < 2; j++) {
                for (int k{0}; k < 2; k++) {
                    vec3 weight_v{u - i, v - j, w - k};
                    accum += (i * uu + (1 - i) * (1 - uu)) *
                             (j * vv + (1 - j) * (1 - vv)) *
                             (k * ww + (1 - k) * (1 - ww)) *
                             dot(c[i][j][k], weight_v);
                }
            }
        }
        return accum;
    }
};
