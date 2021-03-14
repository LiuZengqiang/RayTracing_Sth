#include "initProject.h"
#include <iostream>
#include "globalFunction.h"

#define TINYOBJLOADER_IMPLEMENTATION

#include "Scene.h"

int main() {
    int height = 500;
    int width = 500;
    Scene scene("D:\\CLionProjects\\RayTracing_Sth\\input\\cornellbox\\", "cornellbox.obj");
    return 0;
}
