#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char **argv) {

    // Change the definition here to change resolution
    // Make sure to work with lower resolutions until you are ready to test higher resolutions
    // Scene scene(512, 512);
    Scene scene(256, 256, 2000000);

    Material *red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material *green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material *white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);

    Material *mirror = new Material(MIRROR, Vector3f(0.0f));
    mirror->Kd =  Vector3f(0.14f, 0.45f, 0.091f); //Vector3f(0.95f, 0.93f, 0.88f);
    mirror->ior = 10.f;

    Material *light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f + 0.058f, 0.747f + 0.258f, 0.747f) +
                                             15.6f * Vector3f(0.740f + 0.287f, 0.740f + 0.160f, 0.740f) +
                                             18.4f * Vector3f(0.737f + 0.642f, 0.737f + 0.159f, 0.737f)));
    light->Kd = Vector3f(0.65f);

    // Material* mglassBall = new Material(Transparent);
    // mglassBall->ior_d = 1.5f;
    // mglassBall->SetSmoothness(.9f);
    Material* glass = new Material(TR, Vector3f(0.0f));
    glass->ior = 1.5f;
    glass->roughness = 0.01;
    // glass->Ks = Vector3f(1.f);
    // glass->Kd = Vector3f(0.2f, 0.2f, 0.05f);


    MeshTriangle floor("../models/cornellbox/floor.obj", white);
    MeshTriangle shortbox("../models/cornellbox/shortbox.obj", white);
    MeshTriangle tallbox("../models/cornellbox/tallbox.obj", white);
    MeshTriangle left("../models/cornellbox/left.obj", red);
    MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/light.obj", light);
    //MeshTriangle teapot("../models/cornellbox/teapot.obj", mirror);

    Sphere glassBallLeft(Vector3f(380.0, 95.0, 250.0), 80.0f, glass);
    Sphere glassBallRight(Vector3f(180.0, 95.0, 200.0), 80.0f, glass);

    scene.Add(&floor);
    // scene.Add(&shortbox);
    //scene.Add(&teapot);
    // scene.Add(&tallbox);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);
    scene.Add(&glassBallLeft);
    scene.Add(&glassBallRight);

    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count()
              << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
              << " seconds\n";

    return 0;
}