#pragma once

#include "hittable.h"

#include <memory>
#include <vector>

class hittable_list : public hittable {
  private:
    aabb bbox;

  public:
    std::vector<std::shared_ptr<hittable>> objects{};

  public:
    hittable_list() = default;
    hittable_list(std::shared_ptr<hittable> object) { add(object); }

    void clear() { objects.clear(); }

    void add(std::shared_ptr<hittable> object) {
        objects.push_back(object);
        bbox = aabb(bbox, object->bounding_box());
    }

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override;

    aabb bounding_box() const override { return bbox; }
};

inline bool hittable_list::hit(const ray &r, interval ray_t,
                               hit_record &rec) const {
    hit_record temp_rec;
    hit_record result;
    bool hit_anything{false};
    auto closest_so_far{ray_t.max};

    for (const auto &object : objects) {
        if (object->hit(r, interval{ray_t.min, closest_so_far}, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            result = temp_rec;
        }
    }
    if (hit_anything) {
        rec = result;
    }
    return hit_anything;
}
