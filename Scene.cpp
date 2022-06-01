//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include <cmath>

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const {
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const {
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()) {
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()) {
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum) {
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object *> &objects,
        float &tNear, uint32_t &index, Object **hitObject) {
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }

    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray) const {
    
    // std::cout<<"inside castray"<<std::endl;
    
    Intersection its = intersect(ray);

    if (!its.happened) return Vector3f(0.f);
    
    Vector3f N = normalize(its.normal);
    Vector3f wo = normalize(-ray.direction);

    /*
        TODO: 
        Implement path tracing here.
    */

    if(its.emit.norm() > 0)
        return its.m->getEmission();

    // DIRECT ILLUMINATION - contribution from light source
    // sample light from light source
    Intersection area_p;
    float pdf_light_area;
    sampleLight(area_p, pdf_light_area);

    // convert to solid angle pdf_omega = ||x' - x||^2/Acos(θ')
    float light_obj_dist = ((area_p.coords - its.coords).norm()) * ((area_p.coords - its.coords).norm());
    float area_theta = dotProduct(normalize(its.coords - area_p.coords), normalize(area_p.normal));
    float pdf_light_solidAngle = light_obj_dist * pdf_light_area / area_theta;
    
    // cast ray from its to light
    Vector3f wi_dir = normalize(area_p.coords - its.coords);
    Ray out_dir = Ray(its.coords, wi_dir);
    Intersection hitPoint_dir = intersect(out_dir);
    
    Vector3f L_i = area_p.emit;
    Vector3f f_r_dir = its.m->eval(wi_dir, wo, N);
    float theta_dir = dotProduct(wi_dir, N);

    Vector3f wi = area_p.coords - its.coords;
    // L_dir = L_i * f_r * cos θ / pdf_light(ω from p to x’)
    Vector3f L_dir = Vector3f(0.f);
    if(theta_dir >= cos(M_PI/2) && hitPoint_dir.emit.norm() == area_p.emit.norm())
    //if(theta_dir >= cos(M_PI/2) && wi.norm() < hitPoint_dir.distance)
    //&& dotProduct(-wi_dir, normalize(area_p.normal)) > 0
        L_dir = L_i * f_r_dir * theta_dir / pdf_light_solidAngle;

    // INDIRECT ILLUMINATION - contribution from other reflectors
    Vector3f L_indir = Vector3f(0.f);

    // Russian Roulette
    float ksi = get_random_float();
    if(ksi > RussianRoulette)
        return L_dir + L_indir;

    // Randomly sample the hemisphere toward ωi
    Vector3f wi_indir = normalize(its.m->sample(wo, N)); 

    // cast ray from its to sampled direction
    Ray out_indir = Ray(its.coords, wi_indir);
    Intersection hitPoint_indir = intersect(out_indir);

    float pdf_brdf = its.m->pdf(wi_indir, wo, N);
    Vector3f f_r_indir = its.m->eval(wi_indir, wo, N);
    float theta_indir = dotProduct(wi_indir, N);
    
    if(hitPoint_indir.emit.norm() != area_p.emit.norm()) 
    //if (hitPoint_indir.happened && !hitPoint_indir.obj->hasEmit())
        L_indir = castRay(out_indir) * f_r_indir * theta_indir / pdf_brdf / RussianRoulette;

    return L_dir + L_indir;
}