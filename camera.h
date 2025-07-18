#pragma once

#include <cmath>

#include "color.h"
#include "hittable.h"
#include "material.h"
#include "utils.h"
#include "vec3.h"

class camera {
  public:
    double aspect_ratio{1.0};
    int image_width{100};
    int samples_per_pixel{10};
    int max_depth{10};
    color background;

    double vfov{90};
    point3 lookfrom{0, 0, 0};
    point3 lookat{0, 0, -1};
    vec3 vup{0, 1, 0};

    double defocus_angle{0};
    double focus_dist{10};

  private:
    int image_height;
    double pixel_samples_scale;
    point3 center;
    vec3 pixel_delta_u;
    vec3 pixel_delta_v;
    point3 pixel00_loc;
    vec3 u, v, w;
    vec3 defocus_disk_u;
    vec3 defocus_disk_v;

  public:
    void render(const hittable &world) {
        initialize();
        std::vector buffer(image_height, std::vector<vec3>(image_width));
        for (int j{0}; j < image_height; j++) {

            std::clog << "\rScanlines remaining: " << (image_height - j) << ' '
                      << std::flush;
            // clang-format off
            #pragma omp parallel for schedule(dynamic)
            // clang-format on
            for (int i = 0; i < image_width; i++) {
                color pixel_color{0, 0, 0};
                for (int sample{0}; sample < samples_per_pixel; sample++) {
                    ray r = get_ray(i, j);
                    pixel_color += ray_color(r, max_depth, world);
                }
                buffer[j][i] = pixel_samples_scale * pixel_color;
            }
        }
        std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j{0}; j < image_height; j++) {
            for (int i{0}; i < image_width; i++) {
                const auto &pixel_color{buffer[j][i]};
                write_color(std::cout, pixel_color);
            }
        }
        std::clog << "\rOutputing.            \n";
        std::clog << "\rDone.                 \n";
    }

  private:
    void initialize() {
        image_height =
            static_cast<int>(std::lround(image_width / aspect_ratio));
        image_height = std::max(image_height, 1);

        pixel_samples_scale = 1.0 / samples_per_pixel;

        center = lookfrom;

        auto theta{degrees_to_radians(vfov)};
        auto h{std::tan(theta / 2)};
        auto viewport_height{2.0 * h * focus_dist};
        auto viewport_width{viewport_height *

                            (static_cast<double>(image_width) / image_height)};

        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        vec3 viewport_u{viewport_width * u};
        vec3 viewport_v{viewport_height * -v};

        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        auto viewport_upper_left{center - (focus_dist * w) - viewport_u / 2 -
                                 viewport_v / 2};
        pixel00_loc =
            viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        auto defocus_radius{focus_dist *
                            std::tan(degrees_to_radians(defocus_angle / 2))};
        defocus_disk_u = defocus_radius * u;
        defocus_disk_v = defocus_radius * v;
    }

    ray get_ray(int i, int j) const {
        auto offset{sample_square()};
        auto pixel_sample{pixel00_loc + (i + offset.x()) * pixel_delta_u +
                          (j + offset.y()) * pixel_delta_v};
        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction{pixel_sample - ray_origin};
        auto ray_time = random_double();
        return {ray_origin, ray_direction, ray_time};
    }

    vec3 sample_square() const {
        double random_x{random_double() - 0.5};
        double random_y{random_double() - 0.5};
        return {random_x, random_y, 0};
    }

    color ray_color(const ray &r, int depth, const hittable &world) const {
        if (depth <= 0) {
            return color{0, 0, 0};
        }
        hit_record rec;

        // If the ray hits nothing, return the background color.
        if (!world.hit(r, interval(0.001, infinity), rec))
            return background;

        ray scattered;
        color attenuation;
        color color_from_emission = rec.mat->emitted(rec.u, rec.v, rec.p);

        if (!rec.mat->scatter(r, rec, attenuation, scattered))
            return color_from_emission;

        color color_from_scatter =
            attenuation * ray_color(scattered, depth - 1, world);

        return color_from_emission + color_from_scatter;
    }
    point3 defocus_disk_sample() const {
        auto p{random_in_unit_disk()};
        return center + (p.x() * defocus_disk_u) + (p.y() * defocus_disk_v);
    }
};
