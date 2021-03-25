//
// Created by Sth on 2021/3/24.
//

#ifndef RAYTRACING_STH_TEXTURE_HPP
#define RAYTRACING_STH_TEXTURE_HPP

#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include <string>

class Texture {
public:
    Texture(std::string texture_file_path) : texture_file_path_(texture_file_path) {
        initTexture();
    };

    glm::vec3 getPixel(float u, float v) {
        if (data_ == nullptr) {
            debug::cerrStr("data is nullptr", "Texture::getPixel::");
            return {0.0f, 0.0f, 0.0f};
        }
        if (u > 1.0) {
            u = u - static_cast<int >(u);
        } else if (u < 0.0f) {
            u *= -1;
            u = u - static_cast<int >(u);
            u = 1.0f - u;
        }
        if (v > 1.0) {
            v = v - static_cast<int >(v);
        } else if (v < 0.0f) {
            v *= -1;
            v = v - static_cast<int >(v);
            v = 1.0f - v;
        }


        int x = u * width_;
        int y = v * height_;

        return glm::vec3(static_cast<float >(data_[(y * width_ + x) * 3 + 0]) / 255.0f,
                         static_cast<float >(data_[(y * width_ + x) * 3 + 1]) / 255.0f,
                         static_cast<float >(data_[(y * width_ + x) * 3 + 2]) / 255.0f);
    }

    ~Texture() {
        stbi_image_free(data_);
    }

private:
    void initTexture() {
        data_ = stbi_load(texture_file_path_.c_str(), &width_, &height_, &channel_, 0);
        if (data_ == nullptr) {
            debug::cerrStr(texture_file_path_, "Texture::initTexture::open ", " failed.");
        }
    };

    int width_ = -1;
    int height_ = -1;
    int channel_ = 0;
    std::string texture_file_path_ = "";
    unsigned char *data_ = nullptr;
};

#endif //RAYTRACING_STH_TEXTURE_HPP
