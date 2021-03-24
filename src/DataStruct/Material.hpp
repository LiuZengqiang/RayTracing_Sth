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

    glm::vec3 Kd_ = {0.0f, 0.0f, 0.0f}; // 漫反射

    glm::vec3 Ks_ = {0.0f, 0.0f, 0.0f}; // 镜面反射
    float Ns_ = 0.0f;   // 镜面反射指数
    glm::vec3 Le_ = {0.0f, 0.0f, 0.0f}; // 光

    float Ni_ = 1.0f;   //光密度，折射率。1.0f表示光不会弯曲

    Texture *texture_ = nullptr;

};
//
//inline Ray Material::getWoRay(const glm::vec3 &wi, const glm::vec3 normal, const glm::vec3 ori) {
//
//    if (glm::dot(wi, normal) < 0.0f) {
//        debug::cerrStr("wi and normal has different direction.", "Material::getWoRay");
//    }
//    Ray ret_ray;
//    glm::vec3 wo(0.0f);
//
//    // refraction
//    if (Ni_ > 1.0 + EPSLION) {
//        ret_ray.type = RAY_TYPE::TRANSMISSION;
//        std::cout << "???transmission???" << std::endl;
//
//    } else { //specular and diffuse
//        ret_ray.type = RAY_TYPE::SPECULAR;
//        glm::vec3 reflect_dir = global::getReflection(wi, normal);
//
//        float theta = std::acos(pow(global::getUniform(), 1.0f / (Ns_ + 1)));
//        float phi = global::getUniform() * 2.0f * Pi;
//        float r = std::sin(theta);
//        wo.x = r * std::cos(phi);
//        wo.y = cos(theta);
//        wo.z = r * std::sin(phi);
//        wo = global::toWorld(reflect_dir, wo);
//    }
////    } else {//diffuse
////        // diffuse use normal
////        ret_ray.type = RAY_TYPE::DIFFUSE;
//////        float y = global::clamp(global::getUniform(), 0.0f, 1.0f);
//////        float r = global::clamp(std::sqrt(1.0 - y * y), 0.0f, 1.0f);
//////        float phi = global::getUniform() * 2.0f * Pi;
//////        float x = std::cos(phi) * r;
//////        float z = std::sin(phi) * r;
//////        wo = glm::vec3(x, y, z);
//////        wo = global::toWorld(normal, wo);
////
////        float theta = std::asin(std::sqrt(global::getUniform()));
////        float phi = global::getUniform() * 2.0f * Pi;
////        float r = std::sin(theta);
////        wo.x = r * std::cos(phi);
////        wo.y = cos(theta);
////        wo.z = r * std::sin(phi);
////        wo = global::toWorld(normal, wo);
////    }
//
//    ret_ray.direction = wo;
//    ret_ray.origin = ori;
//    ret_ray.direction_inv = glm::vec3(1.0f / wo.x, 1.0f / wo.y, 1.0f / wo.z);
//    return ret_ray;
//
//    // if has refraction
//    if (this->Ni_ > 1.0 + EPSLION * 2) {
//        // TODO::add refraction
//        ret_ray.type = RAY_TYPE::TRANSMISSION;
//        return ret_ray;
//    } else {    // has not refraction
//
//        float kd_length = glm::length(Kd_);
//        float ks_length = glm::length(Ks_);
//
//        float theta = 0.0f; // theta: 与Y轴的夹角
//
//        glm::vec3 n(0.0f, 1.0f, 0.0f);
//        // 镜面反射
//        if (ks_length > EPSLION && (kd_length / (kd_length + ks_length)) < global::getUniform()) {
//            theta = std::acos(pow(global::getUniform(), 1.0f / (Ns_ + 1)));
//            n = global::getReflection(wi, normal);
//            ret_ray.type = RAY_TYPE::SPECULAR;
//        } else {    // 漫反射 XY平面上cos分布
//            theta = std::asin(std::sqrt(global::getUniform()));
//            ret_ray.type = RAY_TYPE::DIFFUSE;
//        }
//        // phi: XZ平面上与x正方形的夹角，顺时针为正
//        float phi = global::getUniform() * 2.0f * Pi;
//        float r = std::sin(theta);
//        wo.x = r * std::cos(phi);
//        wo.y = cos(theta);
//        wo.z = r * std::sin(phi);
//
//        wo = global::toWorld(n, wo);
//
//        ret_ray.direction = wo;
//        ret_ray.origin = ori;
//        ret_ray.direction_inv = glm::vec3(1.0f / wo.x, 1.0f / wo.y, 1.0f / wo.z);
//        return ret_ray;
//    }
////
////    glm::vec3 wo(0.0f);
////
////    if (this->type_ == MATERIAL_TYPE::DIFFUSE) {
////        // uniform sample on the hemisphere
////        float y = global::dealOutError(global::getUniform(), 0.0f, 1.0f);
////
////        float r = global::dealOutError(std::sqrt(1.0 - y * y), 0.0f, 1.0f);
////
////        float theta = global::getUniform() * 2.0f * Pi;
////        float x = std::cos(theta) * r;
////        float z = std::sin(theta) * r;
////        wo = glm::vec3(x, y, z);
////    } else {
////
////    }
////    // to world
//}
// todo:: add lambertian and glossy material
// lambertian: Ks is 0
// glossy: Ks is not 0
// light: Le is not 0

#endif //RAYTRACING_STH_MATERIAL_HPP
