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
        initScene();
        buildBVHTree();
        bvhTree_->debug();
    };
private:
    void initScene();

    void buildBVHTree();

    std::vector<Triangle *> triangles_;
    std::vector<Material *> materials_;

    std::string file_path_;
    std::string file_name_;
    BVHTree *bvhTree_ = nullptr;
};

#endif //RAYTRACING_STH_SCENE_H
