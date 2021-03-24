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
            debug::cerrStr("data is nullptr", "Texture::getPixel");
            return {0.0f, 0.0f, 0.0f};
        } else if (u < 0 || v < 0 || u > 1.0f || v > 1.0f) {
            debug::cerrStr("out of range", "Texture::getPixel");
            return {0.0f, 0.0f, 0.0f};
        }
        int x = v * width_;
        int y = u * height_;
        return glm::vec3(static_cast<float >(data_[(y * width_ + x) * 3 + 0]),
                         static_cast<float >(data_[(y * width_ + x) * 3 + 1]),
                         static_cast<float >(data_[(y * width_ + x) * 3 + 2]));
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
