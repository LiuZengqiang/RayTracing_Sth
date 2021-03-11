//
// Created by Sth on 2021/2/23.
//

#ifndef RAYTRACING_STH_GLOBALFUNCTION_H
#define RAYTRACING_STH_GLOBALFUNCTION_H

#include "glm/glm.hpp"
#include <vector>
#include <math.h>
#include <iostream>

#define COUT true
#define Pi 3.141592654
#define EPSLION 1e-6
namespace debug {
    inline void coutInt(int i, std::string pre = "", std::string suf = "") {
        if (COUT) {
            std::cout << pre << i << suf << std::endl;
        }
    }

    inline void coutFloat(float i, std::string pre = "", std::string suf = "") {
        if (COUT) {
            std::cout << pre << i << suf << std::endl;
        }
    }

    inline void coutVec3(const glm::vec3 &v, std::string pre = "", std::string suf = "") {
        if (COUT) {
            std::cout << pre << "(" << v.x << "," << v.y << "," << v.z << ")" << suf << std::endl;
        }
    }

    inline void coutVec3s(const std::vector<glm::vec3> &vs, std::string pre = "", std::string suf = "") {
        if (COUT) {
            std::cout << pre;
            for (int i = 0; i < vs.size(); i++) {
                std::cout << "(" << vs[i].x << "," << vs[i].y << "," << vs[i].z << ")";
            }
            std::cout << suf << std::endl;
        }
    }

    inline void coutStr(const std::string str, std::string pre = "", std::string suf = "") {
        if (COUT) {
            std::cout << pre << str << suf << std::endl;
        }
    }

    inline void cerrStr(const std::string str, std::string pre = "", std::string suf = "") {
        if (COUT) {
            std::cerr << pre << str << suf << std::endl;
        }
    }

}
namespace global {
    inline bool floatEqual(const float &a, const float &b) {
        return abs(a - b) < EPSLION;
    }

    /***
     * If a ray has intersection with a triangle whatever in back of front of the triangle.
     * @param v0
     * @param v1
     * @param v2
     * @param ori
     * @param dir
     * @param t
     * @param u
     * @param v
     * @return bool
     */
    inline bool hasIntersectionRayTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &ori,
                                const glm::vec3 &dir, float &t_near, float &u, float &v) {

        glm::vec3 e1 = v1 - v0;
        glm::vec3 e2 = v2 - v0;
        glm::vec3 p = glm::cross(dir, p);
        float det = glm::dot(e1, p);
        if (det <= 0.0f) {
            return false;
        }

        glm::vec3 t = ori - v0;
        u = glm::dot(t, p);

        if (u < 0 || u > det) {
            return false;
        }
        glm::vec3 q = glm::cross(t, e1);
        v = glm::dot(dir, q);
        if (v < 0 || u + v > det) {
            return false;
        }

        float invDet = 1.0f / det;
        t_near = glm::dot(e2, q) * invDet;
        u *= invDet;
        v *= invDet;
        return true;
    }


    inline bool
    hasIntersectionRayBound(const glm::vec3 &p_min, const glm::vec3 &p_max, const glm::vec3 &ori, const glm::vec3 &dir,
                    const glm::vec3 &div_inv) {
        float t_min_x = (p_min.x - ori.x) * div_inv.x;
        float t_max_x = (p_max.x - ori.x) * div_inv.x;
        if (dir.x < 0.0f) {
            std::swap(t_min_x, t_max_x);
        }
        if (t_max_x < 0.0) {
            return false;
        }
        float t_min_y = (p_min.y - ori.y) * div_inv.y;
        float t_max_y = (p_max.y - ori.y) * div_inv.y;
        if (dir.y < 0.0f) {
            std::swap(t_min_y, t_max_y);
        }
        if (t_max_y < 0.0f) {
            return false;
        }
        float t_min_z = (p_min.z - ori.z) * div_inv.z;
        float t_max_z = (p_max.z - ori.z) * div_inv.z;
        if (dir.z < 0.0f) {
            std::swap(t_min_z, t_max_z);
        }
        if (t_max_z < 0.0f) {
            return false;
        }

        float t_enter = std::max(t_min_x, std::max(t_min_y, t_min_z));
        float t_exit = std::min(t_max_x, std::min(t_max_y, t_max_z));
        return t_enter < t_exit;
    };
}

#endif //RAYTRACING_STH_GLOBALFUNCTION_H
