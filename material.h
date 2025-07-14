#pragma once

#include "color.h"
#include "hittable.h"
#include "texture.h"
#include "vec3.h"

class material {
  public:
    virtual ~material() = default;

    virtual bool scatter([[maybe_unused]] const ray &r_in,
                         [[maybe_unused]] const hit_record &rec,
                         [[maybe_unused]] color &attenuation,
                         [[maybe_unused]] ray &scattered) const {
        return false;
    }

    virtual color emitted([[maybe_unused]] double u, [[maybe_unused]] double v,
                          [[maybe_unused]] const point3 &p) const {
        return color(0, 0, 0);
    }
};

class lambertian : public material {
  private:
    std::shared_ptr<texture> tex;

  public:
    lambertian(const color &albedo)
        : tex{std::make_shared<solid_color>(albedo)} {}
    lambertian(std::shared_ptr<texture> tex) : tex{std::move(tex)} {}

    bool scatter([[maybe_unused]] const ray &r_in, const hit_record &rec,
                 color &attenuation, ray &scattered) const override {
        auto scatter_direction = rec.normal + random_unit_vector();

        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;

        scattered = ray{rec.p, scatter_direction, r_in.time()};
        attenuation = tex->value(rec.u, rec.v, rec.p);
        return true;
    }
};

class metal : public material {
  private:
    color albedo;
    double fuzz;

  public:
    metal(const color &albedo, double fuzz)
        : albedo{albedo}, fuzz{fuzz < 1.0 ? fuzz : 1.0} {}

    bool scatter(const ray &r_in, const hit_record &rec, color &attenuation,
                 ray &scattered) const override {
        vec3 reflected = reflect(r_in.direction(), rec.normal);
        reflected = unit_vector(reflected) + fuzz * random_unit_vector();
        scattered = ray{rec.p, reflected, r_in.time()};
        attenuation = albedo;
        return dot(scattered.direction(), rec.normal) > 0;
    }
};

class dielectric : public material {
  private:
    // Refractive index in vacuum or air, or the ratio of the material's
    // refractive index over the refractive index of the enclosing media
    double refraction_index;

  public:
    dielectric(double refraction_index) : refraction_index{refraction_index} {}

    bool scatter(const ray &r_in, const hit_record &rec, color &attenuation,
                 ray &scattered) const override {
        attenuation = color{1.0, 1.0, 1.0};
        double ri =
            rec.front_face ? (1.0 / refraction_index) : refraction_index;

        vec3 unit_direction{unit_vector(r_in.direction())};
        double cos_theta{std::fmin(-dot(unit_direction, rec.normal), 1.0)};
        double sin_theta{std::sqrt(1.0 - cos_theta * cos_theta)};
        bool cannot_refract{ri * sin_theta > 1.0};
        vec3 direction;
        if (cannot_refract || reflectance(cos_theta, ri) > random_double()) {
            direction = reflect(unit_direction, rec.normal);
        } else {
            direction = refract(unit_direction, rec.normal, ri);
        }

        scattered = ray{rec.p, direction, r_in.time()};
        return true;
    }

  private:
    static double reflectance(double cosine, double refraction_index) {
        // Use Schlick's approximation for reflectance.
        auto r0{(1 - refraction_index) / (1 + refraction_index)};
        r0 = r0 * r0;
        return r0 + (1 - r0) * std::pow((1 - cosine), 5);
    }
};

class diffuse_light : public material {
  public:
    diffuse_light(std::shared_ptr<texture> tex) : tex(tex) {}
    diffuse_light(const color &emit)
        : tex(std::make_shared<solid_color>(emit)) {}

    color emitted(double u, double v, const point3 &p) const override {
        return tex->value(u, v, p);
    }

  private:
    std::shared_ptr<texture> tex;
};

class isotropic : public material {
  public:
    isotropic(const color &albedo)
        : tex(std::make_shared<solid_color>(albedo)) {}
    isotropic(std::shared_ptr<texture> tex) : tex(tex) {}

    bool scatter(const ray &r_in, const hit_record &rec, color &attenuation,
                 ray &scattered) const override {
        scattered = ray(rec.p, random_unit_vector(), r_in.time());
        attenuation = tex->value(rec.u, rec.v, rec.p);
        return true;
    }

  private:
    std::shared_ptr<texture> tex;
};

class DisneyMaterial : public material {
  private:
    color base_color;
    double subsurface;
    double metallic;
    double specular;
    double specular_tint;
    double roughness;
    double anisotropic;
    double sheen;
    double sheen_tint;
    double clearcoat;
    double clearcoat_roughness;
    double ior;
    double transmission;
    double transmission_roughness;

  public:
    DisneyMaterial(const color &base_color, double subsurface, double metallic,
                   double specular, double specular_tint, double roughness,
                   double anisotropic, double sheen, double sheen_tint,
                   double clearcoat, double clearcoat_roughness, double ior,
                   double transmission, double transmission_roughness)
        : base_color(base_color), subsurface(subsurface), metallic(metallic),
          specular(specular), specular_tint(specular_tint),
          roughness(roughness), anisotropic(anisotropic), sheen(sheen),
          sheen_tint(sheen_tint), clearcoat(clearcoat),
          clearcoat_roughness(clearcoat_roughness), ior(ior),
          transmission(transmission),
          transmission_roughness(transmission_roughness) {
        // Validate parameters
        roughness = clamp(roughness, 0.001, 1.0);
        clearcoat_roughness = clamp(clearcoat_roughness, 0.001, 1.0);
        transmission_roughness = clamp(transmission_roughness, 0.0, 1.0);
    }

    bool scatter(const ray &r_in, const hit_record &rec, color &attenuation,
                 ray &scattered) const override {
        // Handle transmission
        if (transmission > 0.0 && random_double() < transmission) {
            double refraction_ratio = rec.front_face ? (1.0 / ior) : ior;
            vec3 unit_direction = unit_vector(r_in.direction());

            double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;

            if (cannot_refract ||
                reflectance(cos_theta, refraction_ratio) > random_double()) {
                direction = reflect(unit_direction, rec.normal);
            } else {
                if (transmission_roughness > 0.0) {
                    double alpha =
                        transmission_roughness * transmission_roughness;
                    vec3 h = sample_ggx(rec.normal, alpha, 0.0);
                    direction = refract(unit_direction, h, refraction_ratio);
                    direction = unit_vector(direction);
                } else {
                    direction =
                        refract(unit_direction, rec.normal, refraction_ratio);
                }
            }

            scattered = ray(rec.p, direction, r_in.time());
            attenuation = base_color;
            return true;
        }

        // Regular BSDF evaluation
        vec3 wo = unit_vector(-r_in.direction());
        vec3 wi;
        double pdf;
        color f = evaluate(wo, wi, rec, pdf);

        scattered = ray(rec.p, wi, r_in.time());
        attenuation = f;
        return pdf > 0;
    }

    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }

    static vec3 refract(const vec3 &uv, const vec3 &n, double etai_over_etat) {
        auto cos_theta = fmin(dot(-uv, n), 1.0);
        vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
        vec3 r_out_parallel =
            -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
        return r_out_perp + r_out_parallel;
    }

  private:
    // Helper function to transform from local to world space
    vec3 local_to_world(const vec3 &v, const vec3 &n) const {
        vec3 t = cross(n, vec3(0, 1, 0));
        if (t.length_squared() < 0.001)
            t = cross(n, vec3(1, 0, 0));
        t = unit_vector(t);
        vec3 b = cross(n, t);
        return v.x() * t + v.y() * b + v.z() * n;
    }

    // Sample GGX distribution
    vec3 sample_ggx(const vec3 &n, double alpha, double anisotropic) const {
        if (anisotropic == 0.0) {
            // Isotropic GGX sampling
            double r1 = random_double();
            double r2 = random_double();
            double theta = atan(alpha * sqrt(r1) / sqrt(1.0 - r1));
            double phi = 2.0 * pi * r2;

            vec3 h =
                vec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
            return local_to_world(h, n);
        } else {
            // Anisotropic sampling
            double aspect = sqrt(1.0 - 0.9 * anisotropic);
            double alpha_x = alpha / aspect;
            double alpha_y = alpha * aspect;

            double r1 = random_double();
            double r2 = random_double();

            double phi =
                atan(alpha_y / alpha_x * tan(2.0 * pi * r2 + 0.5 * pi));
            if (r2 > 0.5)
                phi += pi;

            double sin_phi = sin(phi);
            double cos_phi = cos(phi);
            double alpha2 = alpha_x * alpha_x * cos_phi * cos_phi +
                            alpha_y * alpha_y * sin_phi * sin_phi;
            double tan_theta = sqrt(r1 / (1.0 - r1)) * sqrt(alpha2);
            double cos_theta = 1.0 / sqrt(1.0 + tan_theta * tan_theta);
            double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            vec3 h = vec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
            return local_to_world(h, n);
        }
    }

    // Schlick Fresnel approximation
    color fresnel_schlick(double cos_theta, const color &f0) const {
        if (metallic > 0.0) {
            return f0 + (base_color - f0) * std::pow(1.0 - cos_theta, 5.0);
        }
        return f0 +
               (color(1.0, 1.0, 1.0) - f0) * std::pow(1.0 - cos_theta, 5.0);
    }

    // Smith masking-shadowing function
    double smith_g1(const vec3 &v, const vec3 &h, double alpha) const {
        double cos_theta_v = dot(v, h);
        double a = 1.0 / (alpha * std::sqrt(1.0 - cos_theta_v * cos_theta_v));
        return a < 1.6 ? (3.535 * a + 2.181 * a * a) /
                             (1.0 + 2.276 * a + 2.577 * a * a)
                       : 1.0;
    }

    // GGX distribution (with optional anisotropy)
    double ggx_distribution(const vec3 &n, const vec3 &h, double alpha) const {
        double cos_theta = dot(n, h);
        if (cos_theta <= 0)
            return 0.0;

        if (anisotropic == 0.0) {
            // Isotropic case
            double alpha2 = alpha * alpha;
            double denom = cos_theta * cos_theta * (alpha2 - 1.0) + 1.0;
            return alpha2 / (pi * denom * denom);
        } else {
            // Anisotropic case
            vec3 t = cross(n, vec3(0, 1, 0));
            if (t.length_squared() < 0.001)
                t = cross(n, vec3(1, 0, 0));
            t = unit_vector(t);
            vec3 b = cross(n, t);

            double aspect = sqrt(1.0 - 0.9 * anisotropic);
            double alpha_x = alpha / aspect;
            double alpha_y = alpha * aspect;

            double h_dot_t = dot(h, t);
            double h_dot_b = dot(h, b);
            double h_dot_n = dot(h, n);

            double term = (h_dot_t * h_dot_t) / (alpha_x * alpha_x) +
                          (h_dot_b * h_dot_b) / (alpha_y * alpha_y) +
                          h_dot_n * h_dot_n;
            return 1.0 / (pi * alpha_x * alpha_y * term * term);
        }
    }

    // Evaluate Disney BSDF
    color evaluate(const vec3 &wo, vec3 &wi, const hit_record &rec,
                   double &pdf) const {
        // Sample based on material properties
        if (random_double() < 0.5) {
            // Sample diffuse
            wi = random_cosine_direction(rec.normal);
            pdf = 0.5 * abs(dot(wi, rec.normal)) / pi;
        } else {
            // Sample specular
            double alpha = roughness * roughness;
            vec3 h = sample_ggx(rec.normal, alpha, anisotropic);
            wi = reflect(-wo, h);
            if (dot(wi, rec.normal) <= 0)
                return color(0, 0, 0);

            double d = ggx_distribution(rec.normal, h, roughness);
            pdf = 0.5 * d * abs(dot(rec.normal, h)) / (4.0 * abs(dot(wo, h)));
        }

        double n_dot_wi = dot(rec.normal, wi);
        double n_dot_wo = dot(rec.normal, wo);

        if (n_dot_wi <= 0 || n_dot_wo <= 0)
            return color(0, 0, 0);

        // Calculate all components
        vec3 h = unit_vector(wi + wo);

        // Diffuse component
        double energy_bias = mix(0.0, 0.5, subsurface);
        double energy_factor = mix(1.0, 1.0 / 1.51, subsurface);
        double fd90 = energy_bias + 2.0 * roughness * n_dot_wo * n_dot_wo;
        double fi = (1.0 + (fd90 - 1.0) * std::pow(1.0 - n_dot_wi, 5.0));
        double fo = (1.0 + (fd90 - 1.0) * std::pow(1.0 - n_dot_wo, 5.0));
        color diffuse = base_color * fi * fo * n_dot_wi * energy_factor / pi;
        color subsurface =
            base_color * (1.25 / pi) *
            (fi * fo * (1.0 / (n_dot_wi + n_dot_wo) - 0.5) + 0.5);
        color diffuse_comp = mix(diffuse, subsurface, this->subsurface);

        // Specular component
        double dielectric_f0_val = 0.08 * specular * (1.0 - specular_tint);
        color dielectric_f0 =
            color(dielectric_f0_val, dielectric_f0_val, dielectric_f0_val) +
            base_color * specular_tint;
        color f0 = mix(dielectric_f0, base_color, metallic);
        color fresnel = fresnel_schlick(abs(dot(h, wo)), f0);
        double d = ggx_distribution(rec.normal, h, roughness);
        double g = smith_g1(wo, h, roughness) * smith_g1(wi, h, roughness);
        color specular_comp = fresnel * d * g / (4.0 * n_dot_wi * n_dot_wo);

        // Sheen component (velvet-like highlight)
        color sheen_color = mix(color(1.0, 1.0, 1.0), base_color, sheen_tint);
        double sheen_scale = (1.0 - metallic) * sheen;
        color sheen_comp =
            sheen_color * sheen_scale * std::pow(1.0 - abs(dot(h, wo)), 5.0);

        // Clearcoat layer (extra glossy coat)
        double clearcoat_alpha = mix(0.1, 0.001, clearcoat_roughness);
        double clearcoat_d = ggx_distribution(rec.normal, h, clearcoat_alpha);
        double clearcoat_g =
            smith_g1(wo, h, clearcoat_alpha) * smith_g1(wi, h, clearcoat_alpha);
        color clearcoat_fresnel =
            fresnel_schlick(abs(dot(h, wo)), color(0.04, 0.04, 0.04));
        color clearcoat_comp = clearcoat * 0.25 * clearcoat_fresnel *
                               clearcoat_d * clearcoat_g /
                               (4.0 * n_dot_wi * n_dot_wo);

        // Combine all components with energy conservation
        color result =
            (1.0 - metallic) * diffuse_comp * (1.0 - specular) + specular_comp;
        result += sheen_comp * (1.0 - metallic);
        result += clearcoat_comp;
        return result;
    }

    // Linear interpolation between two doubles
    static double mix(double a, double b, double t) {
        return (1.0 - t) * a + t * b;
    }

    // Linear interpolation between two colors
    color mix(const color &a, const color &b, double t) const {
        return (1.0 - t) * a + t * b;
    }

    // Clamp value between min and max
    static double clamp(double val, double min, double max) {
        if (val < min)
            return min;
        if (val > max)
            return max;
        return val;
    }
    vec3 random_cosine_direction(const vec3 &n) const {
        auto r1 = random_double();
        auto r2 = random_double();

        // Sample hemisphere with cosine-weighted distribution
        double phi = 2 * pi * r1;
        double x = cos(phi) * sqrt(r2);
        double y = sin(phi) * sqrt(r2);
        double z = sqrt(1 - r2);

        // Transform from local to world coordinates
        vec3 t, b;
        if (fabs(n.x()) > fabs(n.y())) {
            t = unit_vector(vec3(n.z(), 0, -n.x()));
        } else {
            t = unit_vector(vec3(0, -n.z(), n.y()));
        }
        b = cross(n, t);

        return x * t + y * b + z * n;
    }
};
