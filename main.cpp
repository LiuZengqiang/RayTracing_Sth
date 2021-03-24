#include "initProject.h"
#include <iostream>
#include "globalFunction.h"

#define TINYOBJLOADER_IMPLEMENTATION

#include "Scene.hpp"

int main() {
    int height = 500;
    int width = 500;
    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\cornellbox\\", "cornellbox.obj");
//    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\diningroom\\", "diningroom.obj");
//    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\", "test.obj");
    return 0;
}
