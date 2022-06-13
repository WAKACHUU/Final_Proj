//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"


inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00001;
//const float EPSILON = 0.0001;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // PHOTON MAPPING
    scene.photonMapping();


    // for submission, render 4 images with 1, 4, 16, and 64 spp
    // change the spp value to change number of path samples per pixel
    int spp = 1;
    std::cout << "SPP: " << spp << "\n";
    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5 * rand() / RAND_MAX) / (float)scene.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5 * rand() / RAND_MAX) / (float)scene.height) * scale;

            Vector3f lookat_pos = Vector3f(278 - x, 273 + y, -799);
            Vector3f dir = normalize(lookat_pos - eye_pos);

            for (int k = 0; k < spp; k++){
                framebuffer[m] += scene.castRay(Ray(eye_pos, dir)) / spp;  
            }

            m++;
        }
        UpdateProgress(j / (float)scene.height);
    }
    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        
        float r_level = clamp(0.f, 1.f, pow(framebuffer[i].x, 1.f/2.2f));
        float g_level = clamp(0.f, 1.f, pow(framebuffer[i].y, 1.f/2.2f));
        float b_level = clamp(0.f, 1.f, pow(framebuffer[i].z, 1.f/2.2f));
        color[0] = (unsigned char)(255 * r_level);
        color[1] = (unsigned char)(255 * g_level);
        color[2] = (unsigned char)(255 * b_level);

        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
