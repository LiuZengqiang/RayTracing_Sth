//
// Created by Sth on 2021/3/2.
//

#ifndef RAYTRACING_STH_TRIANGLE_HPP
#define RAYTRACING_STH_TRIANGLE_HPP

#include "glm/glm.hpp"
#include "Material.hpp"
#include "Intersection.hpp"
#include "globalFunction.h"
#include "Ray.hpp"
#include <tuple>
#include <vector>
#include <time.h>

class Bound {
public:
    glm::vec3 p_min_ = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
    glm::vec3 p_max_ = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    glm::vec3 center_ = glm::vec3{0.0f, 0.0f, 0.0f};

    Bound() {
        p_min_ = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
        p_max_ = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
        center_ = glm::vec3{0.0f, 0.0f, 0.0f};
    }

    Bound(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2) {
        setPmin(v0);
        setPmin(v1);
        setPmin(v2);
        setPmax(v0);
        setPmax(v1);
        setPmax(v2);
        center_ = glm::vec3{(p_min_.x + p_min_.x) * 0.5f, (p_min_.y + p_max_.y) * 0.5f, (p_min_.z + p_max_.z) * 0.5f};
    }

    Bound(const glm::vec3 &p1, const glm::vec3 &p2) {
        p_min_ = glm::vec3(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
        p_max_ = glm::vec3(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmin(p1.z, p2.z));
        center_ = glm::vec3{(p_min_.x + p_min_.x) * 0.5f, (p_min_.y + p_max_.y) * 0.5f, (p_min_.z + p_max_.z) * 0.5f};
    }

    Bound(std::vector<glm::vec3> &&vs);

    void setPmin(glm::vec3 p) {
        p_min_.x = std::fmin(p_min_.x, p.x);
        p_min_.y = std::fmin(p_min_.y, p.y);
        p_min_.z = std::fmin(p_min_.z, p.z);
    }

    void setPmax(glm::vec3 p) {
        p_max_.x = std::fmax(p_max_.x, p.x);
        p_max_.y = std::fmax(p_max_.y, p.y);
        p_max_.z = std::fmax(p_max_.z, p.z);
    }

    void debug();

    Bound operator+(Bound b2) const;

    Bound operator+(glm::vec3 p) const;
};

inline Bound Bound::operator+(Bound b2) const {
    Bound ret;
    ret.p_min_ = glm::vec3(fmin(p_min_.x, b2.p_min_.x), fmin(p_min_.y, b2.p_min_.y), fmin(p_min_.z, b2.p_min_.z));
    ret.p_max_ = glm::vec3(fmax(p_max_.x, b2.p_max_.x), fmax(p_max_.y, b2.p_max_.y), fmax(p_max_.z, b2.p_max_.z));
    ret.center_ = glm::vec3{(ret.p_min_.x + ret.p_min_.x) * 0.5f, (ret.p_min_.y + ret.p_max_.y) * 0.5f,
                            (ret.p_min_.z + ret.p_max_.z) * 0.5f};
    return ret;
}

inline Bound Bound::operator+(glm::vec3 p) const {
    Bound ret;
    ret.p_min_ = glm::vec3(fmin(p_min_.x, p.x), fmin(p_min_.y, p.y), fmin(p_min_.z, p.z));
    ret.p_max_ = glm::vec3(fmax(p_max_.x, p.x), fmax(p_max_.y, p.y), fmax(p_max_.z, p.z));
    ret.center_ = glm::vec3{(ret.p_min_.x + ret.p_min_.x) * 0.5f, (ret.p_min_.y + ret.p_max_.y) * 0.5f,
                            (ret.p_min_.z + ret.p_max_.z) * 0.5f};
    return ret;
}

inline Bound::Bound(std::vector<glm::vec3> &&vs) {
    for (size_t i = 0; i < vs.size(); i++) {
        setPmin(vs[i]);
        setPmax(vs[i]);
    }
    center_ = glm::vec3{(p_min_.x + p_min_.x) * 0.5f, (p_min_.y + p_max_.y) * 0.5f, (p_min_.z + p_max_.z) * 0.5f};
}

inline void Bound::debug() {
    debug::coutVec3(p_min_, "p_min:");
    debug::coutVec3(p_max_, "p_max:");
}

class Triangle {
public:
    // construct function
    Triangle(const std::vector<glm::vec3> &ver,
             const std::vector<glm::vec3> &tex,
             const std::vector<glm::vec3> &nor,
             Material *material = nullptr)
            : v0_(ver[0]), v1_(ver[1]), v2_(ver[2]), t0_(tex[0]), t1_(tex[1]), t2_(tex[2]), material_(material) {
        e0_ = v1_ - v0_;
        e1_ = v2_ - v0_;
        if (glm::length(nor[0]) <= EPSLION && glm::length(nor[1]) <= EPSLION && glm::length(nor[2])) {
            normal_ = glm::normalize(glm::cross(e0_, e1_));
        } else {
            normal_ = (nor[0] + nor[1] + nor[2]) / glm::vec3(3.0f, 3.0f, 3.0f);
        }
        area_ = glm::length(glm::dot(e0_, e1_)) * 0.5f;
        bound_ = Bound(v0_, v1_, v2_);
    };

    /***
     * Get a intersection of a ray with the triangle.Must has different direction.
     * @param ray
     * @return
     */
    Intersection getIntersectionWithLimit(const Ray &ray);

    Intersection getIntersectionWithoutLimit(const Ray &ray);

    glm::vec3 getEmission();

    Bound getBound();

    float getArea();

    /***
     * If the triangle is a light, sample a point in the triangle.
     * Using (u,v) coordinate, random sample a value(ux, vx) then the result coordinate is (ux*e0, vx*e1)
     * @return Intersection
     */
    Intersection getSample();

    Material *getMaterial();

    /***
     * If the triangle has emmision
     * @return
     */
    bool hasEmission();

    glm::vec3 v0_, v1_, v2_;
    glm::vec3 t0_, t1_, t2_;
    glm::vec3 e0_, e1_;
    glm::vec3 normal_;
    float area_ = 0.0f;
    Material *material_ = nullptr;
    Bound bound_;
};

inline Bound Triangle::getBound() {
    return bound_;
}

inline float Triangle::getArea() {
    return area_;
}

inline Material *Triangle::getMaterial() {
    return material_;
}

/**
 * get a sample intersection
 * @return Intersection without happened_ information
 */
inline Intersection Triangle::getSample() {
    assert(this->hasEmission());
    Intersection ret_inter;
    float u = global::getUniform();
    float v = global::getUniform();

    if (u + v > 1.0f) {
        float temp_v = v;
        v = 1.0f - u;
        u = 1.0f - temp_v;
    }

    float t = 1.0f - u - v;
    ret_inter.coords_ = glm::vec3(v0_ + e0_ * u + e1_ * v);
    ret_inter.texture_coords_ = glm::vec3(t * t0_ + u * t1_ + v * t2_);
    ret_inter.normal_ = normal_;
    ret_inter.triangle_ = this;
    return ret_inter;
}

inline bool Triangle::hasEmission() {
    return this->material_->hasEmission();
}

inline Intersection Triangle::getIntersectionWithLimit(const Ray &ray) {
    Intersection ret_inter;

    // ray has intersection in triangle's back
    if (glm::dot(ray.direction, normal_) > 0) {
        return ret_inter;
    }
    float u, v, t_near = 0.0f;
    if (global::hasIntersectionRayTriangle(v0_, v1_, v2_, ray.origin, ray.direction, t_near, u, v)) {
        ret_inter.happened_ = true;
        ret_inter.coords_ = ray(t_near);
        ret_inter.texture_coords_ = glm::vec3((1.0f - u - v) * t0_ + u * t1_ + v * t2_);;
        ret_inter.normal_ = normal_;
        ret_inter.distance_ = t_near;
        ret_inter.triangle_ = this;
        return ret_inter;
    } else {
        return ret_inter;
    }
}

inline Intersection Triangle::getIntersectionWithoutLimit(const Ray &ray) {
    Intersection ret_inter;
    float u, v, t_near = 0.0f;
    if (global::hasIntersectionRayTriangle(v0_, v1_, v2_, ray.origin, ray.direction, t_near, u, v)) {
        ret_inter.happened_ = true;
        ret_inter.coords_ = ray(t_near);
        ret_inter.texture_coords_ = glm::vec3((1.0f - u - v) * t0_ + u * t1_ + v * t2_);;
        ret_inter.normal_ = normal_;
        ret_inter.distance_ = t_near;
        ret_inter.triangle_ = this;
        return ret_inter;
    } else {
        return ret_inter;
    }
}

inline glm::vec3 Triangle::getEmission() {
    return this->material_->getEmission();
}

#endif //RAYTRACING_STH_TRIANGLE_HPP
