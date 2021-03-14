//
// Created by Sth on 2021/2/27.
//

#ifndef RAYTRACING_STH_MATERIAL_HPP
#define RAYTRACING_STH_MATERIAL_HPP
// material
#include "globalFunction.h"
#include "glm/glm.hpp"
#include <vector>

enum class MATERIAL_TYPE {
    DIFFUSE,
    LIGHT
};

class Material {
public:
    Material(glm::vec3 Kd = {0.0f, 0.0f, 0.0f},
             glm::vec3 Ks = {0.0f, 0.0f, 0.0f},
             float Ns = 0.0f,
             glm::vec3 Le = {0.0f, 0.0f, 0.0f}
    ) : Kd_(Kd), Ks_(Ks), Ns_(Ns), Le_(Le) {
        if (glm::length(Le_) >= EPSLION) {
            type_ = MATERIAL_TYPE::LIGHT;
        }
        Kd_.y = Kd_.y < EPSLION ? Kd_.x : Kd_.y;
        Kd_.z = Kd_.z < EPSLION ? Kd_.x : Kd_.z;
    };

    // sample a ray by wi and normal
    glm::vec3 getWo(const glm::vec3 &wi, const glm::vec3 normal);

    bool hasEmission();

    glm::vec3 getEmission();

    float getPdf(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &normal);

    glm::vec3 getFr(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &normal);

    void debug() {
        debug::coutVec3(Kd_, "Kd:", "");
        debug::coutVec3(Ks_, "Ks:", "");
        debug::coutFloat(Ns_, "Ns:", "");
        debug::coutVec3(Le_, "Le:", "");
    }

private:

    //Kd 0 1 0
    //Ks 0 0 0
    //Ns 1
    glm::vec3 Kd_ = {0.0f, 0.0f, 0.0f};
    glm::vec3 Ks_ = {0.0f, 0.0f, 0.0f};
    float Ns_ = 0.0f;
    glm::vec3 Le_ = {0.0f, 0.0f, 0.0f};;
    MATERIAL_TYPE type_ = MATERIAL_TYPE::DIFFUSE;
};

inline bool Material::hasEmission() {
    return type_ == MATERIAL_TYPE::LIGHT;
}

inline glm::vec3 Material::getEmission() {
    return Le_;
}

inline glm::vec3 Material::getWo(const glm::vec3 &wi, const glm::vec3 normal) {
    glm::vec3 wo(0.0f);

    if (this->type_ == MATERIAL_TYPE::DIFFUSE) {
        // uniform sample on the hemisphere
        float y = global::dealOutError(global::getUniform(), 0.0f, 1.0f);

        float r = global::dealOutError(std::sqrt(1.0 - y * y), 0.0f, 1.0f);

        float theta = global::getUniform() * 2 * Pi;
        float x = std::cos(theta) * r;
        float z = std::sin(theta) * r;
        wo = glm::vec3(x, y, z);
    } else {

    }
    // to world
    return global::toWorld(normal, wo);
}

inline float Material::getPdf(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &normal) {
    switch (this->type_) {
        case MATERIAL_TYPE::DIFFUSE:
            return glm::dot(wo, normal) > 0.0f ? 0.5f / Pi : 0.0f;
    }
    return 0;
}

inline glm::vec3 Material::getFr(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &normal) {
    switch (this->type_) {
        case MATERIAL_TYPE::DIFFUSE:
            return glm::dot(wo, normal) > 0.0f ? this->Kd_ * glm::vec3(1.0f / Pi) : glm::vec3(0.0f);
            break;
    }
    debug::cerrStr("Error:");
    return glm::vec3(0.0f);
}

#endif //RAYTRACING_STH_MATERIAL_HPP
