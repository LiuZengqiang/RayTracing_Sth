//
// Created by Sth on 2021/2/27.
//

#ifndef RAYTRACING_STH_MATERIAL_HPP
#define RAYTRACING_STH_MATERIAL_HPP
// material
#include "globalFunction.h"
#include "glm/glm.hpp"
#include "Ray.hpp"
#include "Texture.hpp"
#include <vector>

class Material {
public:
    Material() {};

    Material(glm::vec3 Kd = {0.0f, 0.0f, 0.0f},
             glm::vec3 Ks = {0.0f, 0.0f, 0.0f},
             float Ns = 0.0f,
             glm::vec3 Le = {0.0f, 0.0f, 0.0f},
             std::string name = "none",
             Texture *texture = nullptr
    ) : Kd_(Kd), Ks_(Ks), Ns_(Ns), Le_(Le), material_name_(name), texture_(texture) {
    };
    // sample a ray by wi and normal
    // Ray getWoRay(const glm::vec3 &wi, const glm::vec3 normal, const glm::vec3 ori);

    Ray getWoRayDiffuse(const glm::vec3 &wi, const glm::vec3 normal, const glm::vec3 ori) {
        Ray ret_ray;
//        ret_ray.type = RAY_TYPE::DIFFUSE;
        glm::vec3 wo(0.0f);
        float u = global::getUniform();
        float theta = std::acos(1.0 - u);
        float phi = global::getUniform() * 2.0f * Pi;
        float r = std::sin(theta);
        wo.x = r * std::cos(phi);
        wo.y = std::cos(theta);
        wo.z = r * std::sin(phi);
        wo = global::toWorld(normal, wo);
        ret_ray.origin = ori;
        ret_ray.direction = wo;
        ret_ray.direction_inv = 1.0f / ret_ray.direction;
        return ret_ray;
    };

    Ray getWoRaySpecular(const glm::vec3 &wi, const glm::vec3 normal, const glm::vec3 ori) {
        Ray ret_ray;
        glm::vec3 wo(0.0f);
//        ret_ray.type = RAY_TYPE::SPECULAR;
        glm::vec3 reflect_dir = global::getReflection(wi, normal);
        float u = global::getUniform();
        float theta = std::acos(pow(u, 1.0f / (Ns_ + 1)));
        float phi = global::getUniform() * 2.0f * Pi;
        float r = std::sin(theta);
        wo.x = r * std::cos(phi);
        wo.y = cos(theta);
        wo.z = r * std::sin(phi);
        wo = global::toWorld(reflect_dir, wo);
        ret_ray.origin = ori;
        ret_ray.direction = wo;
        ret_ray.direction_inv = 1.0f / ret_ray.direction;
        return ret_ray;
    };

    bool hasEmission() {
        return glm::length(Le_) > EPSLION;
    };

    glm::vec3 getEmission() {
        return Le_;
    };

    glm::vec3 getKd() {
        return Kd_;
    };

    glm::vec3 getKs() {
        return Ks_;
    };

    float getNs() {
        return Ns_;
    }

    Texture *getTexture() {
        return texture_;
    }

    std::string getName() {
        return material_name_;
    };

    void debug() {
        debug::coutVec3(Kd_, "Kd:", "");
        debug::coutVec3(Ks_, "Ks:", "");
        debug::coutFloat(Ns_, "Ns:", "");
        debug::coutVec3(Le_, "Le:", "");
    }

private:
    std::string material_name_ = "None";

    glm::vec3 Kd_ = {0.0f, 0.0f, 0.0f}; // diffuse

    glm::vec3 Ks_ = {0.0f, 0.0f, 0.0f}; // specular
    float Ns_ = 0.0f;   // 镜面反射指数
    glm::vec3 Le_ = {0.0f, 0.0f, 0.0f}; // emission

//    glm::vec3 Ka_ = {0.5f, 0.5f, 0.5f}; // ambient

    float Ni_ = 1.0f;   //光密度，折射率。1.0f表示光不会弯曲

    Texture *texture_ = nullptr;

};

#endif //RAYTRACING_STH_MATERIAL_HPP
