#include <cassert>
#include <memory>

#include "bvh.h"
#include "camera.h"
#include "constant_medium.h"
#include "hittable_list.h"
#include "material.h"
#include "quad.h"
#include "sphere.h"
#include "texture.h"

void bouncing_spheres() {
    hittable_list world;

    auto checker = std::make_shared<checker_texture>(0.32, color(.2, .3, .1),
                                                     color(.9, .9, .9));

    world.add(std::make_shared<sphere>(point3(0, -1000, 0), 1000,
                                       std::make_shared<lambertian>(checker)));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2,
                          b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                std::shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = std::make_shared<lambertian>(albedo);
                    auto center2 = center + vec3{0, random_double(0, 0.5), 0};
                    world.add(std::make_shared<sphere>(center, center2, 0.2,
                                                       sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = std::make_shared<metal>(albedo, fuzz);
                    world.add(
                        std::make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = std::make_shared<dielectric>(1.5);
                    world.add(
                        std::make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = std::make_shared<dielectric>(1.5);
    world.add(std::make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(std::make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(std::make_shared<sphere>(point3(4, 1, 0), 1.0, material3));
    world = hittable_list(std::make_shared<bvh_node>(world));

    camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;
    cam.background = color(0.70, 0.80, 1.00);

    cam.vfov = 20;
    cam.lookfrom = point3{13, 2, 3};
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0.6;
    cam.focus_dist = 10.0;
    cam.render(world);
}

void checkered_spheres() {
    hittable_list world;

    auto checker = std::make_shared<checker_texture>(0.32, color(.2, .3, .1),
                                                     color(.9, .9, .9));

    world.add(std::make_shared<sphere>(point3(0, -10, 0), 10,
                                       std::make_shared<lambertian>(checker)));
    world.add(std::make_shared<sphere>(point3(0, 10, 0), 10,
                                       std::make_shared<lambertian>(checker)));

    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;
    cam.background = color(0.70, 0.80, 1.00);

    cam.vfov = 20;
    cam.lookfrom = point3(13, 2, 3);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(world);
}

void earth() {

    auto earth_texture = std::make_shared<image_texture>("earthmap.jpg");
    auto earth_surface = std::make_shared<lambertian>(earth_texture);
    auto globe = std::make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);

    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;
    cam.background = color(0.70, 0.80, 1.00);

    cam.vfov = 20;
    cam.lookfrom = point3(0, 0, 12);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(hittable_list(globe));
}

void perlin_spheres() {
    hittable_list world;

    auto pertext = std::make_shared<noise_texture>(4);
    world.add(std::make_shared<sphere>(point3(0, -1000, 0), 1000,
                                       std::make_shared<lambertian>(pertext)));
    world.add(std::make_shared<sphere>(point3(0, 2, 0), 2,
                                       std::make_shared<lambertian>(pertext)));

    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;
    cam.background = color(0.70, 0.80, 1.00);

    cam.vfov = 20;
    cam.lookfrom = point3(13, 2, 3);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(world);
}

void quads() {
    hittable_list world;

    // Materials
    auto left_red = std::make_shared<lambertian>(color(1.0, 0.2, 0.2));
    auto back_green = std::make_shared<lambertian>(color(0.2, 1.0, 0.2));
    auto right_blue = std::make_shared<lambertian>(color(0.2, 0.2, 1.0));
    auto upper_orange = std::make_shared<lambertian>(color(1.0, 0.5, 0.0));
    auto lower_teal = std::make_shared<lambertian>(color(0.2, 0.8, 0.8));

    // Quads
    world.add(std::make_shared<quad>(point3(-3, -2, 5), vec3(0, 0, -4),
                                     vec3(0, 4, 0), left_red));
    world.add(std::make_shared<quad>(point3(-2, -2, 0), vec3(4, 0, 0),
                                     vec3(0, 4, 0), back_green));
    world.add(std::make_shared<quad>(point3(3, -2, 1), vec3(0, 0, 4),
                                     vec3(0, 4, 0), right_blue));
    world.add(std::make_shared<quad>(point3(-2, 3, 1), vec3(4, 0, 0),
                                     vec3(0, 0, 4), upper_orange));
    world.add(std::make_shared<quad>(point3(-2, -3, 5), vec3(4, 0, 0),
                                     vec3(0, 0, -4), lower_teal));

    camera cam;

    cam.aspect_ratio = 1.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;
    cam.background = color(0.70, 0.80, 1.00);

    cam.vfov = 80;
    cam.lookfrom = point3(0, 0, 9);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(world);
}
void simple_light() {
    hittable_list world;

    auto pertext = std::make_shared<noise_texture>(4);
    world.add(std::make_shared<sphere>(point3(0, -1000, 0), 1000,
                                       make_shared<lambertian>(pertext)));
    world.add(std::make_shared<sphere>(point3(0, 2, 0), 2,
                                       make_shared<lambertian>(pertext)));

    auto difflight = std::make_shared<diffuse_light>(color(4, 4, 4));
    world.add(std::make_shared<sphere>(point3(0, 7, 0), 2, difflight));
    world.add(std::make_shared<quad>(point3(3, 1, -2), vec3(2, 0, 0),
                                     vec3(0, 2, 0), difflight));

    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;
    cam.background = color(0, 0, 0);

    cam.vfov = 20;
    cam.lookfrom = point3(26, 3, 6);
    cam.lookat = point3(0, 2, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(world);
}
void cornell_box() {
    hittable_list world;

    auto red = std::make_shared<lambertian>(color(.65, .05, .05));
    auto white = std::make_shared<lambertian>(color(.73, .73, .73));
    auto green = std::make_shared<lambertian>(color(.12, .45, .15));
    auto light = std::make_shared<diffuse_light>(color(15, 15, 15));

    world.add(std::make_shared<quad>(point3(555, 0, 0), vec3(0, 555, 0),
                                     vec3(0, 0, 555), green));
    world.add(std::make_shared<quad>(point3(0, 0, 0), vec3(0, 555, 0),
                                     vec3(0, 0, 555), red));
    world.add(std::make_shared<quad>(point3(343, 554, 332), vec3(-130, 0, 0),
                                     vec3(0, 0, -105), light));
    world.add(std::make_shared<quad>(point3(0, 0, 0), vec3(555, 0, 0),
                                     vec3(0, 0, 555), white));
    world.add(std::make_shared<quad>(point3(555, 555, 555), vec3(-555, 0, 0),
                                     vec3(0, 0, -555), white));
    world.add(std::make_shared<quad>(point3(0, 0, 555), vec3(555, 0, 0),
                                     vec3(0, 555, 0), white));
    std::shared_ptr<hittable> box1 =
        box(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = std::make_shared<rotate_y>(box1, 15);
    box1 = std::make_shared<translate>(box1, vec3(265, 0, 295));
    world.add(box1);

    std::shared_ptr<hittable> box2 =
        box(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = std::make_shared<rotate_y>(box2, -18);
    box2 = std::make_shared<translate>(box2, vec3(130, 0, 65));
    world.add(box2);

    camera cam;

    cam.aspect_ratio = 1.0;
    cam.image_width = 600;
    cam.samples_per_pixel = 200;
    cam.max_depth = 50;
    cam.background = color(0, 0, 0);

    cam.vfov = 40;
    cam.lookfrom = point3(278, 278, -800);
    cam.lookat = point3(278, 278, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(world);
}

void cornell_smoke() {
    hittable_list world;

    auto red = std::make_shared<lambertian>(color(.65, .05, .05));
    auto white = std::make_shared<lambertian>(color(.73, .73, .73));
    auto green = std::make_shared<lambertian>(color(.12, .45, .15));
    auto light = std::make_shared<diffuse_light>(color(7, 7, 7));

    world.add(std::make_shared<quad>(point3(555, 0, 0), vec3(0, 555, 0),
                                     vec3(0, 0, 555), green));
    world.add(std::make_shared<quad>(point3(0, 0, 0), vec3(0, 555, 0),
                                     vec3(0, 0, 555), red));
    world.add(std::make_shared<quad>(point3(113, 554, 127), vec3(330, 0, 0),
                                     vec3(0, 0, 305), light));
    world.add(std::make_shared<quad>(point3(0, 555, 0), vec3(555, 0, 0),
                                     vec3(0, 0, 555), white));
    world.add(std::make_shared<quad>(point3(0, 0, 0), vec3(555, 0, 0),
                                     vec3(0, 0, 555), white));
    world.add(std::make_shared<quad>(point3(0, 0, 555), vec3(555, 0, 0),
                                     vec3(0, 555, 0), white));

    std::shared_ptr<hittable> box1 =
        box(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = std::make_shared<rotate_y>(box1, 15);
    box1 = std::make_shared<translate>(box1, vec3(265, 0, 295));

    std::shared_ptr<hittable> box2 =
        box(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = std::make_shared<rotate_y>(box2, -18);
    box2 = std::make_shared<translate>(box2, vec3(130, 0, 65));

    world.add(std::make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));
    world.add(std::make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));

    camera cam;

    cam.aspect_ratio = 1.0;
    cam.image_width = 600;
    cam.samples_per_pixel = 200;
    cam.max_depth = 50;
    cam.background = color(0, 0, 0);

    cam.vfov = 40;
    cam.lookfrom = point3(278, 278, -800);
    cam.lookat = point3(278, 278, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    cam.render(world);
}

void disney_material_scene() {
    hittable_list world;

    // Ground plane
    auto ground = std::make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(std::make_shared<sphere>(point3(0, -1000, 0), 1000, ground));

    // Create a sphere for each Disney parameter
    const int cols = 7;
    const double spacing = 2.0;
    const double sphere_radius = 0.3;

    // Row 1: Core parameters
    for (int col = 0; col < cols; col++) {
        point3 center((col - cols / 2.0) * spacing, sphere_radius, -spacing);
        color base_color = color(0.8, 0.1, 0.1);
        auto mat = std::make_shared<DisneyMaterial>(
            base_color,
            col == 0 ? 1.0 : 0.0,                    // Subsurface
            col == 1 ? 1.0 : 0.0,                    // Metallic
            col == 2 ? 1.0 : 0.5,                    // Specular
            col == 3 ? 1.0 : 0.0,                    // Specular tint
            col == 4 ? 0.1 : (col == 5 ? 0.9 : 0.5), // Roughness
            col == 6 ? 1.0 : 0.0,                    // Anisotropic
            0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0);
        world.add(std::make_shared<sphere>(center, sphere_radius, mat));
    }

    // Row 2: Additional parameters
    for (int col = 0; col < cols; col++) {
        point3 center((col - cols / 2.0) * spacing, sphere_radius, spacing);
        color base_color = color(0.1, 0.1, 0.8);
        auto mat = std::make_shared<DisneyMaterial>(
            base_color, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0,
            col == 0 ? 1.0 : 0.0,                    // Sheen
            col == 1 ? 1.0 : 0.0,                    // Sheen tint
            col == 2 ? 1.0 : 0.0,                    // Clearcoat
            col == 3 ? 0.1 : (col == 4 ? 0.9 : 0.5), // Clearcoat roughness
            1.5,
            col == 5 ? 1.0 : 0.0,  // Transmission
            col == 6 ? 1.0 : 0.0); // Transmission roughness
        world.add(std::make_shared<sphere>(center, sphere_radius, mat));
    }

    // Add central sphere with all parameters at mid values
    auto default_mat = std::make_shared<DisneyMaterial>(
        color(0.5, 0.5, 0.5), 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        1.5, 0.5, 0.5);
    world.add(std::make_shared<sphere>(point3(0, sphere_radius, 0),
                                       sphere_radius, default_mat));

    // Lights
    auto light = std::make_shared<diffuse_light>(color(15, 15, 15));
    world.add(std::make_shared<quad>(point3(-5, 5, -5), vec3(10, 0, 0),
                                     vec3(0, 0, 10), light));

    camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 800;
    cam.samples_per_pixel = 400;
    cam.max_depth = 50;
    cam.background = color(0.1, 0.1, 0.1);

    cam.vfov = 40;
    cam.lookfrom = point3(0, 5, 10);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;
    cam.render(world);
}

int main() {
    constexpr int CHOICE = 9;
    switch (CHOICE) {
    case 1:
        bouncing_spheres();
        break;
    case 2:
        checkered_spheres();
        break;
    case 3:
        earth();
        break;
    case 4:
        perlin_spheres();
        break;
    case 5:
        quads();
        break;
    case 6:
        simple_light();
        break;
    case 7:
        cornell_box();
        break;
    case 8:
        cornell_smoke();
        break;
    case 9:
        disney_material_scene();
        break;
    default:
        assert(false);
    }
    return 0;
}
