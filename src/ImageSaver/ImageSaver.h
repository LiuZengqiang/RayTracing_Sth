//
// Created by Sth on 2021/2/27.
//

#ifndef RAYTRACING_STH_IMAGESAVER_H
#define RAYTRACING_STH_IMAGESAVER_H

#include "glm/glm.hpp"
#include <vector>
#include <string>

class ImageSaver {
public:
    static ImageSaver *getInstance() {
        if (nullptr == instance_) {
            instance_ = new ImageSaver();
        }
        return instance_;
    }

    /***
     * Save vector data to image(.ppm), the data[0][0] at left up corner of image, data[m][n] at right down corner.
     * @param image_data
     * @param width
     * @param height
     * @param output_name
     */
    void saveImage(std::vector<glm::vec3> &image_data, int width, int height,
                   std::string output_name = "..\\output\\output.ppm") {
        std::cout << "Begin save image." << std::endl;
        FILE *fp = fopen(output_name.c_str(), "wb");
        assert(fp != nullptr);
        (void) fprintf(fp, "P6\n%d %d\n255\n", width, height);

        for (auto i = 0; i < width * height; ++i) {
            static unsigned char color[3];
            color[0] = (unsigned char) (255 * std::pow(global::clamp(image_data[i].x, 0.0f, 1.0f), 0.45f));
            color[1] = (unsigned char) (255 * std::pow(global::clamp(image_data[i].y, 0.0f, 1.0f), 0.45f));
            color[2] = (unsigned char) (255 * std::pow(global::clamp(image_data[i].z, 0.0f, 1.0f), 0.45f));
            fwrite(color, 1, 3, fp);
        }
        fclose(fp);
        std::cout << "End save image." << std::endl;
    }

    static ImageSaver *instance_;
private:
};

ImageSaver *ImageSaver::instance_ = nullptr;
#endif //RAYTRACING_STH_IMAGESAVER_H
