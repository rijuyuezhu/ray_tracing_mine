#pragma once

#include <memory>

#include "hittable.h"
#include "hittable_list.h"
#include "material.h"
#include "vec3.h"

class quad : public hittable {
  private:
    point3 Q;
    vec3 u, v;
    vec3 w;
    std::shared_ptr<material> mat;
    aabb bbox;
    vec3 normal;
    double D;

  public:
    quad(const point3 &Q, const vec3 &u, const vec3 &v,
         std::shared_ptr<material> mat)
        : Q(Q), u(u), v(v), mat(mat) {
        auto n{cross(u, v)};
        normal = unit_vector(n);
        D = dot(normal, Q);
        w = n / dot(n, n);
        set_bounding_box();
    }

    virtual void set_bounding_box() {
        aabb bbox_diagonal1{Q, Q + u + v};
        aabb bbox_diagonal2{Q + u, Q + v};
        bbox = aabb(bbox_diagonal1, bbox_diagonal2);
    }

    aabb bounding_box() const override { return bbox; }

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override {
        auto denom{dot(normal, r.direction())};
        if (std::abs(denom) < 1e-8) {
            return false;
        }
        auto t{(D - dot(normal, r.origin())) / denom};
        if (!ray_t.contains(t)) {
            return false;
        }
        auto intersection{r.at(t)};
        vec3 planar_hitpt_vector = intersection - Q;
        auto alpha = dot(w, cross(planar_hitpt_vector, v));
        auto beta = dot(w, cross(u, planar_hitpt_vector));

        if (!is_interior(alpha, beta, rec))
            return false;
        rec.t = t;
        rec.p = intersection;
        rec.mat = mat;
        rec.set_face_normal(r, normal);

        return true;
    }
    virtual bool is_interior(double a, double b, hit_record &rec) const {
        interval unit_interval = interval(0, 1);
        // Given the hit point in plane coordinates, return false if it is
        // outside the primitive, otherwise set the hit record UV coordinates
        // and return true.

        if (!unit_interval.contains(a) || !unit_interval.contains(b))
            return false;

        rec.u = a;
        rec.v = b;
        return true;
    }
};

inline std::shared_ptr<hittable_list> box(const point3 &a, const point3 &b,
                                          std::shared_ptr<material> mat) {
    // Returns the 3D box (six sides) that contains the two opposite vertices a
    // & b.

    auto sides = std::make_shared<hittable_list>();

    // Construct the two opposite vertices with the minimum and maximum
    // coordinates.
    auto min = point3(std::fmin(a.x(), b.x()), std::fmin(a.y(), b.y()),
                      std::fmin(a.z(), b.z()));
    auto max = point3(std::fmax(a.x(), b.x()), std::fmax(a.y(), b.y()),
                      std::fmax(a.z(), b.z()));

    auto dx = vec3(max.x() - min.x(), 0, 0);
    auto dy = vec3(0, max.y() - min.y(), 0);
    auto dz = vec3(0, 0, max.z() - min.z());

    sides->add(std::make_shared<quad>(point3(min.x(), min.y(), max.z()), dx, dy,
                                      mat)); // front
    sides->add(std::make_shared<quad>(point3(max.x(), min.y(), max.z()), -dz,
                                      dy, mat)); // right
    sides->add(std::make_shared<quad>(point3(max.x(), min.y(), min.z()), -dx,
                                      dy, mat)); // back
    sides->add(std::make_shared<quad>(point3(min.x(), min.y(), min.z()), dz, dy,
                                      mat)); // left
    sides->add(std::make_shared<quad>(point3(min.x(), max.y(), max.z()), dx,
                                      -dz, mat)); // top
    sides->add(std::make_shared<quad>(point3(min.x(), min.y(), min.z()), dx, dz,
                                      mat)); // bottom

    return sides;
}
