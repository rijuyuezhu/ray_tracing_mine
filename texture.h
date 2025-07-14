#pragma once

#include <memory>

#include "color.h"
#include "perlin.h"
#include "rtw_stb_image.h"
#include "vec3.h"

class texture {
  public:
    virtual ~texture() = default;
    virtual color value(double u, double v, const point3 &p) const = 0;
};

class solid_color : public texture {
  private:
    color albedo;

  public:
    solid_color(const color &albedo) : albedo(albedo) {}

    solid_color(double red, double green, double blue)
        : albedo(color(red, green, blue)) {}

    color value([[maybe_unused]] double u, [[maybe_unused]] double v,
                [[maybe_unused]] const point3 &p) const override {
        return albedo;
    }
};

class checker_texture : public texture {
  private:
    double inv_scale;
    std::shared_ptr<texture> even;
    std::shared_ptr<texture> odd;

  public:
    checker_texture(double scale, std::shared_ptr<texture> even,
                    std::shared_ptr<texture> odd)
        : inv_scale(1.0 / scale), even(std::move(even)), odd(std::move(odd)) {}

    checker_texture(double scale, const color &c1, const color &c2)
        : checker_texture(scale, std::make_shared<solid_color>(c1),
                          std::make_shared<solid_color>(c2)) {}

    color value([[maybe_unused]] double u, [[maybe_unused]] double v,
                const point3 &p) const override {
        auto xInterger{static_cast<int>(std::floor(inv_scale * p.x()))};
        auto yInterger{static_cast<int>(std::floor(inv_scale * p.y()))};
        auto zInterger{static_cast<int>(std::floor(inv_scale * p.z()))};

        bool isEven = (xInterger + yInterger + zInterger) % 2 == 0;

        return isEven ? even->value(u, v, p) : odd->value(u, v, p);
    }
};

class image_texture : public texture {
  private:
    rtw_image image;

  public:
    image_texture(const char *filename) : image(filename) {}

    color value(double u, double v,
                [[maybe_unused]] const point3 &p) const override {
        u = interval(0, 1).clamp(u);
        v = 1.0 - interval(0, 1).clamp(v);

        int i{static_cast<int>(u * image.width())};
        int j{static_cast<int>(v * image.height())};
        auto pixel{image.pixel_data(i, j)};
        auto color_scale = 1.0 / 255.0;
        return color(color_scale * pixel[0], color_scale * pixel[1],
                     color_scale * pixel[2]);
    }
};

class noise_texture : public texture {
  private:
    perlin noise;
    double scale;

  public:
    noise_texture(double scale) : scale(scale) {}
    color value([[maybe_unused]] double u, [[maybe_unused]] double v,
                const point3 &p) const override {
        return color(.5, .5, .51) *
               (1 + std::sin(scale * p.z() + 10 * noise.turb(p, 7)));
    }
};
