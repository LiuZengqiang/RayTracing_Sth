//
// Created by Sth on 2021/2/24.
//
// BVH-Tree
// BVH:(Bounding Volume Hierarchy)
// leaf node: list of object, box
// internal node: box, two children
//TODO::Complete BVH

#ifndef RAYTRACING_STH_BVHTREE_HPP
#define RAYTRACING_STH_BVHTREE_HPP

#include "Triangle.hpp"
#include <vector>
#include <glm/glm.hpp>
#include <algorithm>
#include <queue>

enum class SORT_TYPE {
    X, Y, Z
};

struct BVHNode {
    std::vector<Triangle *> triangles_;
    BVHNode *left;
    BVHNode *right;
    Bound bound;
    int triangle_size_ = 0; // just for debug

    BVHNode(std::vector<Triangle *> triangles) {
        triangles_ = triangles;
        left = right = nullptr;
        if (!triangles.empty()) {
            bound = Bound(std::vector<glm::vec3>{triangles[0]->v0_, triangles_[0]->v1_, triangles_[0]->v2_});
            for (size_t i = 1; i < triangles.size(); i++) {
                bound = bound +
                        Bound(std::vector<glm::vec3>{triangles[i]->v0_, triangles_[i]->v1_, triangles_[i]->v2_});
            }
            triangle_size_ = triangles.size();
        }
    }

    bool isLeaf() {
        return left == nullptr && right == nullptr;
    };
};

class BVHTree {
public:
    BVHTree(std::vector<Triangle *> triangles, int max_num_pre_node = 1) : max_num_pre_node_(max_num_pre_node) {

        debug::coutStr("Begin build BVHTree.", "BVHTree::BVHTree()::");
        root_ = buildBVHTree(triangles);
        debug::coutStr("End build BVHTree.", "BVHTree::BVHTree()::");
    }

    BVHNode *getRoot();

    Intersection getIntersection(Ray &ray);

    void debug();

private:
    BVHNode *buildBVHTree(std::vector<Triangle *> triangles_, SORT_TYPE sort_type = SORT_TYPE::X);

    Intersection getIntersection(BVHNode *root, Ray &ray);

    int max_num_pre_node_ = 1;
    BVHNode *root_ = nullptr;
};

inline BVHNode *BVHTree::buildBVHTree(std::vector<Triangle *> triangles, SORT_TYPE sort_type) {
    int triangle_size = triangles.size();
    if (triangles.empty()) {
        return nullptr;
    }

    if (triangles.size() <= max_num_pre_node_) {
        BVHNode *node = new BVHNode(triangles);
        node->triangle_size_ = triangles.size();
        return node;
    } else {
        // x->y->z
        switch (sort_type) {
            case SORT_TYPE::X:
                std::sort(triangles.begin(), triangles.end(),
                          [](const Triangle *a, const Triangle *b) {
                              return (a->bound_.center_.x) < (b->bound_.center_.x);
                          });
                sort_type = SORT_TYPE::Y;
                break;
            case SORT_TYPE::Y:
                std::sort(triangles.begin(), triangles.end(),
                          [](const Triangle *a, const Triangle *b) {
                              return (a->bound_.center_.y) < (b->bound_.center_.y);
                          });
                sort_type = SORT_TYPE::Z;
                break;
            default:
                std::sort(triangles.begin(), triangles.end(),
                          [](const Triangle *a, const Triangle *b) {
                              return (a->bound_.center_.z) < (b->bound_.center_.z);
                          });
                sort_type = SORT_TYPE::X;
        }

        BVHNode *node = new BVHNode(std::vector<Triangle *>());

        std::vector<Triangle *>::iterator mid_iter = triangles.begin() + triangles.size() / 2;

        std::vector<Triangle *> left_triangle(triangles.begin(), mid_iter);
        node->left = buildBVHTree(std::vector<Triangle *>(triangles.begin(), mid_iter), sort_type);

        node->right = buildBVHTree(std::vector<Triangle *>(mid_iter, triangles.end()), sort_type);

        node->bound = node->left->bound + node->right->bound;
        node->triangle_size_ = node->left->triangle_size_ + node->right->triangle_size_;
        return node;
    }
}

inline BVHNode *BVHTree::getRoot() {
    return root_;
}

inline void BVHTree::debug() {
    debug::coutStr("============");
    debug::coutStr("BVH Tree Level traversal.", "BVHTree::");
    std::queue<BVHNode *> que;
    que.push(root_);

    while (!que.empty()) {
        std::queue<BVHNode *> temp_que;
        while (!que.empty()) {
            BVHNode *temp = que.front();
            debug::coutInt(que.front()->triangle_size_);

            debug::coutVec3(temp->bound.p_min_);
            debug::coutVec3(temp->bound.p_max_);
            que.pop();
            if (temp->left) {
                temp_que.push(temp->left);
            }
            if (temp->right) {
                temp_que.push(temp->right);
            }
        }
        std::cout << std::endl;
        que = std::move(temp_que);
    }
    debug::coutStr("============");
}


inline Intersection BVHTree::getIntersection(Ray &ray) {
    return getIntersection(root_, ray);
}

inline Intersection BVHTree::getIntersection(BVHNode *root, Ray &ray) {

    if (root == nullptr ||
        global::hasIntersectionRayBound(root->bound.p_min_, root->bound.p_max_, ray.origin, ray.direction,
                                        ray.direction_inv) == false) {
        return Intersection();
    }
    // should check every triangle in the bvh node
    if (root->isLeaf()) {
        Intersection ret_inter = root->triangles_[0]->getIntersection(ray);
        for (int i = 1; i < root->triangles_.size(); i++) {
            Intersection temp_inter = root->triangles_[i]->getIntersection(ray);
            if (ret_inter.distance_ > temp_inter.distance_) {
                ret_inter = std::move(temp_inter);
            }
        }
        return ret_inter;
    } else {

        Intersection left_inter = getIntersection(root->left, ray);
        Intersection right_inter = getIntersection(root->right, ray);
        return left_inter.distance_ < right_inter.distance_ ? left_inter : right_inter;
    }

}

#endif //RAYTRACING_STH_BVHTREE_HPP
