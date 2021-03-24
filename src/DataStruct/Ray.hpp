//
// Created by Sth on 2021/2/24.
//
// Ray Datastruct
#ifndef RAYTRACING_STH_RAY_HPP
#define RAYTRACING_STH_RAY_HPP

#include <glm/glm.hpp>

enum class RAY_TYPE {
    DIFFUSE,
    SPECULAR,
    TRANSMISSION,
    NONE
};

// Ray struct
struct Ray {
//    enum RAY_TYPE type = RAY_TYPE::NONE;
    glm::vec3 origin = glm::vec3(0.0f);
    glm::vec3 direction = glm::vec3(0.0f);
    glm::vec3 direction_inv;    // for multiply operator is faster  than divide
    Ray(glm::vec3 ori, glm::vec3 dir) : origin(ori), direction(glm::normalize(dir)) {
        direction_inv = glm::vec3(1.0f / direction.x, 1.0f / direction.y, 1.0f / direction.z);
    }

    Ray() {};

    glm::vec3 operator()(double t) const {
        glm::vec3 temp_vec3(direction);
        temp_vec3 *= t;
        return origin + temp_vec3;
    }
};

#endif //RAYTRACING_STH_RAY_HPP
