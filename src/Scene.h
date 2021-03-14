//
// Created by Sth on 2021/2/23.
//

// Material.h .cpp
// Triangle.h .cpp
#ifndef RAYTRACING_STH_SCENE_H
#define RAYTRACING_STH_SCENE_H

#include "tiny_obj_loader.h"
#include "DataStruct/Triangle.hpp"
#include "DataStruct/Material.hpp"
#include "DataStruct/BVHTree.hpp"
#include <vector>
#include <string>

class Scene {
public:
    Scene(std::string file_path, std::string file_name) {
        file_path_ = file_path;
        file_name_ = file_name;
        debug::coutStr("Begin initScene", "Scene::");
        initScene();
        debug::coutStr("Begin buildBVHTree", "Scene::");
        buildBVHTree();
        debug();

        Ray ray;
        pathTracing(ray);
    };

    glm::vec3 pathTracing(const Ray &ray) const;

    void debug();

private:
    glm::vec3 shade(Intersection &p, glm::vec3 wo) const;

    void initScene();

    void buildBVHTree();

    Triangle *sampleLightTriangle() const;

    std::vector<Triangle *> triangles_;
    std::vector<Triangle *> light_triangles_;
    std::vector<Material *> materials_;

    std::string file_path_;
    std::string file_name_;
    BVHTree *bvhTree_ = nullptr;

    float light_area_ = 0.0f;

    float russian_roulette_ = 0.8f;
};

#endif //RAYTRACING_STH_SCENE_H
