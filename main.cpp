#include "initProject.h"
#include <iostream>
#include "globalFunction.h"

#define TINYOBJLOADER_IMPLEMENTATION

#include "Scene.hpp"

int main() {
    int height = 500;
    int width = 500;

//    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\car\\", "car.obj", glm::vec3(8.22f, -0.61f, -9.80f),
//                glm::vec3(7.514f, -0.702f, -9.097f), glm::vec3(-0.065f, 0.996f, 0.065f), 45.0f);
    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\cornellbox\\", "cornellbox.obj");
//    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\diningroom\\", "diningroom.obj");
//    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\", "test.obj");
    return 0;
}
