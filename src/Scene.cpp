//
// Created by Sth on 2021/2/23.
//
#include "globalFunction.h"
#include "Scene.h"
#include <iostream>
#include <string>
#include <sstream>
#include <map>

void Scene::initScene() {
    debug::coutStr("Begin initScene.", "Scene::");

    tinyobj::ObjReaderConfig reader_config;
    reader_config.triangulate = true;

    reader_config.mtl_search_path = file_path_;

    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(file_path_ + file_name_, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "ERROR::Scene::initScene:" << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        debug::coutStr("WARNING::Scene::initScene:" + reader.Warning());
//        std::cout << "WARNING::Scene::initScene:" << reader.Warning();
    }
    auto &attrib = reader.GetAttrib();  // attributes
    auto &shapes = reader.GetShapes();  // shapes
    auto &materials = reader.GetMaterials();    // materials
    debug::coutInt(materials.size(), "material size:", ".");

    // load material
    for (int i = 0; i < materials.size(); i++) {
        // Kd:diffuse
        // Ks:specular
        // Ns:shininess
        // Le:unknown_parameter
        glm::vec3 Kd(0.0f);
        glm::vec3 Ks(0.0f);
        int Ns = 0.0f;
        glm::vec3 Le(0.0f);
        Kd = glm::vec3(materials[i].diffuse[0], materials[i].diffuse[1],
                       materials[i].diffuse[2]);
        Ks = glm::vec3(materials[i].specular[0], materials[i].specular[1],
                       materials[i].specular[2]);
        Ns = float(materials[i].shininess);

        if (materials[i].unknown_parameter.find("Le") !=
            materials[i].unknown_parameter.end()) {
            std::string le = materials[i].unknown_parameter.at("Le");
            std::stringstream ss(le);
            ss >> Le[0] >> Le[1] >> Le[2];
            ss.clear();
        }
        Material *new_material = new Material(Kd, Ks, Ns, Le);
        materials_.push_back(new_material);
    }

    debug::coutStr("shape size" + std::to_string(shapes.size()));
    int face_cnt = 0;
    // loop over shapes
    for (size_t s = 0; s < shapes.size(); ++s) {
        // every shape has a mesh(face)

        // the index of every vertex belong to one mesh in together(shapes[s].mesh.indices[x]).

        size_t index_offset = 0;
        // shapes[s].mesh.num_face_vertices contain the num of vertices per face
        // so the size of .num_face_vertices is the num of face per shape
        // loop over faces
        face_cnt += shapes[s].mesh.num_face_vertices.size();
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); ++f) {

            int v_size = shapes[s].mesh.num_face_vertices[f];
            assert(v_size == 3);
            std::vector<glm::vec3> ver(3, glm::vec3(0.0f));
            std::vector<glm::vec3> tex(3, glm::vec3(0.0f));
            std::vector<glm::vec3> nor(3, glm::vec3(0.0f));
            // loop over vertices
            for (size_t v = 0; v < v_size; ++v) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vertex_x = attrib.vertices[3 * idx.vertex_index + 0];
                tinyobj::real_t vertex_y = attrib.vertices[3 * idx.vertex_index + 1];
                tinyobj::real_t vertex_z = attrib.vertices[3 * idx.vertex_index + 2];
                ver[v] = glm::vec3(vertex_x, vertex_y, vertex_z);

                if (3 * idx.normal_index < attrib.normals.size()) {
                    tinyobj::real_t normal_x = attrib.normals[3 * idx.normal_index + 0];
                    tinyobj::real_t normal_y = attrib.normals[3 * idx.normal_index + 1];
                    tinyobj::real_t normal_z = attrib.normals[3 * idx.normal_index + 2];
                    nor[v] = glm::vec3(normal_x, normal_y, normal_z);
                }

                if (2 * idx.texcoord_index < attrib.texcoords.size()) {
                    tinyobj::real_t texture_x = attrib.texcoords[2 * idx.texcoord_index + 0];
                    tinyobj::real_t texture_y = attrib.texcoords[2 * idx.texcoord_index + 1];
                    tex[v] = glm::vec3(texture_x, texture_y, 0.0f);
                }

            }

            index_offset += v_size;
            // per-face material
            int f_material_idx = shapes[s].mesh.material_ids[f];
            Material *ma_p = (f_material_idx <= -1) ? nullptr : materials_[f_material_idx];

            Triangle *new_triangle = new Triangle(ver, tex, nor, ma_p);
            debug::coutVec3s(ver, "New Triangle");
            triangles_.push_back(new_triangle);
        }
    }

    debug::coutInt(face_cnt, "tiny_load_obj face size:", ".");
    debug::coutInt(triangles_.size(), "triangles face size:", ".");
    debug::coutInt(materials_.size(), "material size:", ".");
}

void Scene::buildBVHTree() {
    bvhTree_ = new BVHTree(triangles_, 1);
    if (bvhTree_->getRoot() == nullptr) {
        debug::cerrStr("root is nullptr!", "buildBVHTree::", "");
    }
    debug::coutInt(bvhTree_->getRoot()->triangle_size_, "root->triangle_size:", "");
}