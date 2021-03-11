//
// Created by Sth on 2021/2/27.
//

#ifndef RAYTRACING_STH_INTERSECTION_H
#define RAYTRACING_STH_INTERSECTION_H
// Intersection
#include "glm/glm.hpp"

class Triangle;

class Intersection {
public:

    bool happened_;
    glm::vec3 coords_;
    glm::vec3 texture_coords_;
    glm::vec3 normal_;
    float distance_;
    Triangle *triangle_;

    Intersection() {
        happened_ = false;
        coords_ = glm::vec3(0.0f, 0.0f, 0.0f);
        texture_coords_ = glm::vec3(0.0f, 0.0f, 0.0f);
        normal_ = glm::vec3(0.0f, 0.0f, 0.0f);
        distance_ = FLT_MAX;
        triangle_ = nullptr;
    }

};

#endif //RAYTRACING_STH_INTERSECTION_H
