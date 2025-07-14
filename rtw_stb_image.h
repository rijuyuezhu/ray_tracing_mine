#pragma once

#define STB_IMAGE_IMPLEMENTATION
#define STBI_FAILURE_USERMSG
#include "external/stb_image.h"

#include <cstdlib>
#include <iostream>
#include <string>

class rtw_image {
  private:
    static constexpr int bytes_per_pixel{3};
    float *fdata{nullptr};
    std::uint8_t *bdata{nullptr};
    int image_width{8};
    int image_height{8};
    int bytes_per_scanline{0};

  public:
    rtw_image() {}
    rtw_image(const char *image_filename) {
        std::string filename{image_filename};
        auto imagedir{std::getenv("RTW_IMAGES")};

        if (imagedir && load(std::string(imagedir) + "/" + filename))
            return;
        if (load(filename))
            return;
        if (load("images/" + filename))
            return;
        if (load("../images/" + filename))
            return;
        if (load("../../images/" + filename))
            return;
        if (load("../../../images/" + filename))
            return;
        if (load("../../../../images/" + filename))
            return;
        if (load("../../../../../images/" + filename))
            return;
        if (load("../../../../../../images/" + filename))
            return;
        std::cerr << "ERROR: Could not load image file '" << image_filename
                  << "'.\n";
        std::exit(EXIT_FAILURE);
    }
    ~rtw_image() {
        delete[] bdata;
        STBI_FREE(fdata);
    }
    bool load(const std::string &filename) {
        auto n{bytes_per_pixel};
        fdata = stbi_loadf(filename.c_str(), &image_width, &image_height, &n,
                           bytes_per_pixel);
        if (fdata == nullptr)
            return false;
        bytes_per_scanline = image_width * bytes_per_pixel;
        convert_to_bytes();
        return true;
    }
    int width() const { return image_width; }
    int height() const { return image_height; }

    const unsigned char *pixel_data(int x, int y) const {
        x = clamp(x, 0, image_width);
        y = clamp(y, 0, image_height);
        return bdata + y * bytes_per_scanline + x * bytes_per_pixel;
    }

  private:
    static int clamp(int x, int low, int high) {
        if (x < low)
            return low;
        else if (x < high)
            return x;
        else
            return high - 1;
    }

    static std::uint8_t float_to_byte(float f) {
        if (f < 0.0f)
            return 0;
        else if (f > 1.0f)
            return 255;
        else
            return static_cast<std::uint8_t>(f * 256.0f);
    }
    void convert_to_bytes() {
        int total_bytes = image_width * image_height * bytes_per_pixel;
        bdata = new std::uint8_t[total_bytes];

        auto bptr{bdata};
        auto fptr{fdata};
        for (auto i{0}; i < total_bytes; i++, bptr++, fptr++) {
            *bptr = float_to_byte(*fptr);
        }
    }
};
