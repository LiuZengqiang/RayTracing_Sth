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

//        debug();
    };

    // sample a ray by wi and normal
    glm::vec3 getWo(const glm::vec3 &wi, const glm::vec3 normal);

    bool hasEmission();

    glm::vec3 getEmission();

    float getPdf(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &normal);

    float getFr(const glm::vec3 &wi, const glm::vec3 &wo, const glm::vec3 &normal);

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

#endif //RAYTRACING_STH_MATERIAL_HPP
