#pragma once

#include <algorithm>

#include "aabb.h"
#include "hittable.h"
#include "hittable_list.h"
#include "utils.h"

class bvh_node : public hittable {
  private:
    std::shared_ptr<hittable> left;
    std::shared_ptr<hittable> right;
    aabb bbox;

  public:
    bvh_node(hittable_list list)
        : bvh_node(list.objects, 0, list.objects.size()) {}

    bvh_node(std::vector<std::shared_ptr<hittable>> &objects, std::size_t start,
             std::size_t end) {
        bbox = aabb::empty;
        for (std::size_t object_index = start; object_index < end;
             object_index++)
            bbox = aabb(bbox, objects[object_index]->bounding_box());
        int axis = bbox.longest_axis();

        auto comparator = [axis](const std::shared_ptr<hittable> &a,
                                 const std::shared_ptr<hittable> &b) {
            return a->bounding_box().axis_interval(axis).min <
                   b->bounding_box().axis_interval(axis).min;
        };

        std::size_t object_span = end - start;
        if (object_span == 1) {
            left = right = objects[start];
        } else if (object_span == 2) {
            left = objects[start];
            right = objects[start + 1];
        } else {
            std::sort(objects.begin() + start, objects.begin() + end,
                      comparator);
            auto mid = start + object_span / 2;
            left = std::make_shared<bvh_node>(objects, start, mid);
            right = std::make_shared<bvh_node>(objects, mid, end);
        }
    }

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override {
        if (!bbox.hit(r, ray_t)) {
            return false;
        }
        bool hit_left = left->hit(r, ray_t, rec);
        bool hit_right = right->hit(
            r, interval(ray_t.min, hit_left ? rec.t : ray_t.max), rec);
        return hit_left || hit_right;
    }

    aabb bounding_box() const override { return bbox; }
};
