//
// Created by Sth on 2021/2/23.
//

// Material.h .cpp
// Triangle.h .cpp
#ifndef RAYTRACING_STH_SCENE_HPP
#define RAYTRACING_STH_SCENE_HPP

#include "globalFunction.h"
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include "ImageSaver/ImageSaver.h"
#include "tiny_obj_loader.h"
#include "DataStruct/Triangle.hpp"
#include "DataStruct/Material.hpp"
#include "DataStruct/BVHTree.hpp"
#include "DataStruct/Texture.hpp"
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
        saveImage();
    };

    void saveImage() {
        image_width_ = 100;
        image_height_ = 100;

        image_data_.resize(image_width_ * image_height_);

        glm::vec3 position(0.0f, 0.0f, 2.5f);
        glm::vec3 lookAt(0.0f, 0.0f, 0.0f);

        glm::vec3 direction = glm::normalize(lookAt - position);
        glm::vec3 up(0.0f, 1.0f, 0.0f);
        up = glm::normalize(glm::cross(glm::cross(direction, up), direction));
        glm::vec3 right = glm::normalize(glm::cross(direction, up));
        float fov = 60.0f;

        fov = fov / 180.0f * Pi;
        float d = image_height_ * 0.5 / std::tan(fov * 0.5);
        glm::vec3 o = position + d * direction;
        glm::vec3 left_up =
                o + 0.5f * image_height_ * up - 0.5f * image_width_ * right;
        glm::vec3 right_up =
                o + 0.5f * image_height_ * up + 0.5f * image_width_ * right;
        glm::vec3 left_bottom =
                o - 0.5f * image_height_ * up - 0.5f * image_width_ * right;

        glm::vec3 du = (right_up - left_up);
        glm::vec3 dv = (left_bottom - left_up);
        debug::coutVec3(left_up, "left_up:");
        debug::coutVec3(right_up, "right_up:");
        debug::coutVec3(left_bottom, "left_bottom:");
//    Ray ray(position, glm::normalize(left_bottom - position));
        Ray ray(position, glm::vec3(0.0f, 0.0f, -1.0f));
        glm::vec3 val = pathTracing(ray);
        debug::coutVec3(val, "val:");

        for (int i = 0; i < image_height_; i++) {
            for (int j = 0; j < image_width_; j++) {


                glm::vec3 val(0, 0, 0);
                for (int k = 0; k < 40; k++) {
                    glm::vec3 target(left_up + ((static_cast<float >(j) + global::getUniform()) / image_width_) * du +
                                     ((static_cast<float >(i) + global::getUniform()) / image_height_) * dv);
                    ray = Ray(position, glm::normalize(target - position));
                    val += pathTracing(ray);
                }
                val /= 40;
                image_data_[i * image_width_ + j] = val;
            }
            debug::coutInt(i);
        }
        ImageSaver *imageSaver = ImageSaver::getInstance();
        imageSaver->saveImage(image_data_, image_width_, image_height_);
    }

//    glm::vec3 pathTracing(const Ray &ray) const;

    glm::vec3 pathTracing(const Ray &ray, int depth = 0) const {
        Intersection intersection = bvhTree_->getIntersection(ray);

        if (intersection.happened_ == false || glm::dot(intersection.normal_, ray.direction) > 0.0f) {
            return glm::vec3(0.0f, 0.0f, 0.0f);
        } else if (intersection.triangle_->hasEmission()) {
            return global::clamp(intersection.triangle_->getEmission());
        } else if (depth <= 10) {
            Material *material = intersection.triangle_->material_;

            glm::vec3 indirect_light(0.0f);

            if (glm::length(material->getKd()) >= EPSLION) {
                Ray next_ray_diffuse = intersection.triangle_->material_->getWoRayDiffuse(-1.0f * ray.direction,
                                                                                          intersection.normal_,
                                                                                          intersection.coords_);
//            next_ray_diffuse.origin = next_ray_diffuse.origin + 0.001f * next_ray_diffuse.direction;
                Intersection next_inter = bvhTree_->getIntersection(next_ray_diffuse);
                if (next_inter.happened_ == false || next_inter.triangle_->hasEmission() ||
                    glm::dot(next_inter.normal_, next_ray_diffuse.direction) > 0.0f) { ;
                } else {
                    if (intersection.triangle_->material_->getTexture() != nullptr) {
                        indirect_light += material->getKd() * pathTracing(next_ray_diffuse, depth + 1) *
                                          material->getTexture()->getPixel(intersection.texture_coords_[0],
                                                                           intersection.texture_coords_[1]);
                    } else {
                        indirect_light += material->getKd() * pathTracing(next_ray_diffuse, depth + 1);
                    }
                }
            }
            if (glm::length(material->getKs()) >= EPSLION) {
                Ray next_ray_specular = intersection.triangle_->material_->getWoRaySpecular(-1.0f * ray.direction,
                                                                                            intersection.normal_,
                                                                                            intersection.coords_);
//            next_ray_specular.origin = next_ray_specular.origin + 0.001f * next_ray_specular.direction;
                Intersection next_inter = bvhTree_->getIntersection(next_ray_specular);
                if (next_inter.happened_ == false || next_inter.triangle_->hasEmission() ||
                    glm::dot(next_inter.normal_, next_ray_specular.direction) > 0.0f) { ;
                } else {
                    indirect_light += material->getKs() * pathTracing(next_ray_specular, depth + 1);
                }
            }
            indirect_light = global::clamp(indirect_light);

//        Ray next_ray = intersection.triangle_->material_->getWoRay(-1.0f * ray.direction, intersection.normal_,
//                                                                   intersection.coords_);
//        // 下一个交点
//        next_ray.origin = next_ray.origin + 0.001f * next_ray.direction;
//
//        Intersection next_inter = bvhTree_->getIntersection(next_ray);
//
//        if (next_inter.happened_ == false || next_inter.triangle_->hasEmission() ||
//            glm::dot(next_inter.normal_, next_ray.direction) > 0.0f) { ;
//        } else {
//            if (next_ray.type != RAY_TYPE::NONE) {
//                indirect_light = pathTracing(next_ray, depth + 1);
//                indirect_light *= (material->getKd() + material->getKs()) * 0.5f;
////                switch (next_ray.type) {
////                    case RAY_TYPE::DIFFUSE:
////                        indirect_light = material->getKd() * indirect_light;
////                        break;
////                    case RAY_TYPE::SPECULAR:
//////                        debug::coutStr("specular", "pathTracing::");
////                        indirect_light = material->getKs() * indirect_light;
////                        break;
////                    case RAY_TYPE::TRANSMISSION:
////                        debug::coutStr("transmission.");
////                        break;
////                    default:
////                        break;
////                }
//            }
////            indirect_light = global::clamp(indirect_light);
//        }
            glm::vec3 direct_light = externalLight(ray, intersection);
//        int sample_size = 10;
//        for (int i = 0; i < light_triangles_.size(); i++) {
//            glm::vec3 single_ret(0.0f);
//
//            for (int j = 0; j < sample_size; j++) {
//                Intersection inter = light_triangles_[i]->getSample();
//
//                Ray direct_ray = Ray(intersection.coords_, glm::normalize(inter.coords_ - intersection.coords_));
//
//                Intersection inter_other = bvhTree_->getIntersection(direct_ray);
//
//                float ray_distance = glm::length(inter.coords_ - intersection.coords_);
//
//                if (inter_other.distance_ < ray_distance ||
//                    glm::dot(intersection.normal_, direct_ray.direction) < 0.0f ||
//                    glm::dot(inter.normal_, direct_ray.direction) > 0.0f) {
//                    continue;
//                } else {
//                    float cos_i = global::clamp(glm::dot(intersection.normal_, direct_ray.direction));
//                    float cos_o = global::clamp(glm::dot(inter.normal_, -1.0f * direct_ray.direction));
//                    glm::vec3 e = cos_i * cos_o / (ray_distance * ray_distance) * light_triangles_[i]->getArea() *
//                                  light_triangles_[i]->getEmission();
//                    single_ret += material->getKd() * e;
//                }
//            }
//
//            single_ret /= sample_size;
//        }
//        std::cout << "direct_light:";
//        debug::coutVec3(direct_light);
            return direct_light + indirect_light;
        } else {
            return glm::vec3(0.0f);
        }
    };

    glm::vec3 externalLight(const Ray &ray, Intersection intersection) const {
        Material *material = intersection.triangle_->material_;
        glm::vec3 ret(0.0f);
        for (int i = 0; i < light_triangles_.size(); i++) {
            glm::vec3 local_rgb(0.0f);
            for (int j = 0; j < 10; j++) {
                Intersection inter = light_triangles_[i]->getSample();
                glm::vec3 light_origin = inter.coords_;
                glm::vec3 light_direction = light_origin - intersection.coords_;
                float length_light = glm::length(light_direction);
                Ray light_ray = Ray(intersection.coords_, glm::normalize(light_direction));

                Intersection inter_other = bvhTree_->getIntersection(light_ray);

                if (inter_other.distance_ >= length_light - EPSLION) {
                    glm::vec3 light_d_n = glm::normalize(light_direction);
                    float cos_i = std::max(0.0f, glm::dot(intersection.normal_, light_d_n));
                    float cos_o = std::max(0.0f, glm::dot(light_triangles_[i]->normal_, -1.0f * light_d_n));
                    float f0 = cos_i * cos_o / (length_light * length_light);
                    glm::vec3 energy = f0 * light_triangles_[i]->getArea() * light_triangles_[i]->getEmission();
                    float kd_length = glm::length(material->getKd());
                    float ks_length = glm::length(material->getKs());

                    if (kd_length > EPSLION) {
                        float cos_i_o = global::clamp(glm::dot(-ray.direction, intersection.normal_));
                        if (material->getTexture() != nullptr) {
                            local_rgb +=
                                    material->getKd() * energy * (cos_i_o) * (kd_length / (kd_length + ks_length)) *
                                    material->getTexture()->getPixel(intersection.texture_coords_[0],
                                                                     intersection.texture_coords_[1]);
                        } else {
                            local_rgb += material->getKd() * energy * (cos_i_o) * (kd_length / (kd_length + ks_length));
                        }

                    }

                    if (ks_length > EPSLION) {
                        glm::vec3 reflect_wi = global::getReflection(light_ray.direction, intersection.normal_);
                        float cos_n = global::clamp(
                                std::pow(glm::dot(reflect_wi, -1.0f * ray.direction), material->getNs()));
                        local_rgb += material->getKs() * energy * cos_n * (ks_length / (kd_length + ks_length));
                    }

//                debug::coutVec3(energy, "energy");
//                if (glm::length(material->getKd()) > EPSLION) {
//                    float cos_i_o = glm::dot(light_d_n, intersection.normal_);
//                    local_rgb += material->getKd() * energy;
//                    if (cos_i_o > EPSLION) {
//                        local_rgb += cos_i_o * material->getKd() * energy;
//                    }
//                }
                }
            }
//        debug::coutVec3(local_rgb, "local_rgb");
            ret += local_rgb / 10.0f / Pi;
        }
        global::clamp(ret);
        return ret;
    };

    void debug() {
        std::cout << "Triangle size:" << triangles_.size() << std::endl;
        std::cout << "Light triangle size:" << light_triangles_.size() << std::endl;
        std::cout << "Material size:" << materials_.size() << std::endl;
    };

private:

    void initScene() {


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
            glm::vec3 Kd(0.0f);
            glm::vec3 Ks(0.0f);
            int Ns = 0.0f;
            glm::vec3 Le(0.0f);
            Kd = glm::vec3(materials[i].diffuse[0], materials[i].diffuse[1],
                           materials[i].diffuse[2]);
            Ks = glm::vec3(materials[i].specular[0], materials[i].specular[1],
                           materials[i].specular[2]);
            Ns = float(materials[i].shininess);
            std::string material_name = materials[i].name;

            if (materials[i].unknown_parameter.find("Le") !=
                materials[i].unknown_parameter.end()) {
                std::string le = materials[i].unknown_parameter.at("Le");
                std::stringstream ss(le);
                debug::coutStr(le, "le:");
                ss >> Le[0] >> Le[1] >> Le[2];
                debug::coutVec3(Le, "le:");
                ss.clear();
            }
            // map name
            //        materials[i].diffuse_texname;
            Texture *new_texture = nullptr;


            if (materials[i].diffuse_texname != "") {
                new_texture = new Texture(file_path_ + materials[i].diffuse_texname);
                texture_data_.push_back(new_texture);
            }

            Material *new_material = new Material(Kd, Ks, Ns, Le, material_name, new_texture);
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
                triangles_.push_back(new_triangle);
                if (new_triangle->material_->hasEmission()) {
                    light_triangles_.push_back(new_triangle);
                    light_area_ += new_triangle->getArea();
                }
            }
        }
    };

    void buildBVHTree() {
        bvhTree_ = new BVHTree(triangles_, 1);
        if (bvhTree_->getRoot() == nullptr) {
            debug::cerrStr("root is nullptr!", "buildBVHTree::", "");
        }
    };

    Triangle *sampleLightTriangle() const;

    std::vector<Triangle *> triangles_;
    std::vector<Triangle *> light_triangles_;
    std::vector<Material *> materials_;

    std::string file_path_;
    std::string file_name_;
    BVHTree *bvhTree_ = nullptr;

    float light_area_ = 0.0f;

    float russian_roulette_ = 0.8f;

    int image_width_ = 500;
    int image_height_ = 500;
    std::vector<glm::vec3> image_data_;

    std::vector<Texture *> texture_data_;
};

#endif //RAYTRACING_STH_SCENE_HPP
