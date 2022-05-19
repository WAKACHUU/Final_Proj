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
    Intersection its = intersect(ray);
    if (!its.happened) return Vector3f(0.f);
    
    Vector3f N = normalize(its.normal);
    Vector3f wo = normalize(-ray.direction);

    /*
        TODO: 
        Implement path tracing here.
    */
    if (its.obj->hasEmit()) {
        return its.emit;
    }

    Intersection light;
    float pdf_a;
    sampleLight(light, pdf_a); //now we have light intersect and its pdf
    Vector3f wi = light.coords - its.coords; // wi
    Vector3f wi_n = normalize(wi);
    Vector3f n_light = normalize(light.normal); //N of the light
    Vector3f wi_rev = -wi_n;

    float cosine_light = dotProduct(wi_rev, n_light);
    
    float light_area = 1.0 / pdf_a;
    float pdf_w = wi.norm()*wi.norm() / (cosine_light * light_area);


    Vector3f f_r = its.m->eval(wi_n, wo, N); //BRDF of light

    float cosine = dotProduct(wi_n, N);

    Vector3f L_dir = Vector3f(), L_indir = Vector3f();

    //if it is blocked
    Ray ray_i = Ray(its.coords, wi_n);
    Intersection shdIts = intersect(ray_i); 

    if (!shdIts.happened && cosine_light > 0) { 
        L_dir = light.emit * f_r * cosine / pdf_w; 
        
    }
    else {
        /*wi.norm() < shdIts.distance*/
        if ( wi.norm() >= shdIts.distance || cosine_light <= 0 || dotProduct(wi_n, N) < 0) {
            L_dir = Vector3f(0);
            //std::cout << "!!" << cosine_light << std::endl;
        }
        else if (cosine_light>0) {
            L_dir = light.emit * f_r * cosine / pdf_w; 
        }
    }   
    

//indirect part!
    Vector3f wi_b = normalize(its.m->sample(wo, N));
    float pdf_brdf = its.m->pdf(wi_b, wo, N);

    Ray newIts = Ray(its.coords, wi_b);

    Vector3f f_r_new = its.m->eval(wi_b, wo, N);

    float cosine_new = dotProduct(wi_b, N);

    if (intersect(newIts).happened && !intersect(newIts).obj->hasEmit()) { //not emitting 
        L_indir = castRay(newIts) * f_r_new * cosine_new / pdf_brdf / RussianRoulette;
    }

    return L_dir + L_indir;
}