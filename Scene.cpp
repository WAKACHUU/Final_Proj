//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include <cmath>
#include <typeinfo>

Scene::Scene(int w, int h, int numPhotons) {
    width = w;
    height = h;
    this->numPhotons = numPhotons;
    causticsMap = new PhotonMap(numPhotons);
    globalMap = new PhotonMap(numPhotons);
}

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

void Scene::tracePhoton(Ray &ray, Vector3f &power) const{
    
    // std::cout<<"inside trace photons"<<std::endl;

    // currently only implementing diffuse case

    Intersection hit = intersect(ray);

    // photon doesn't land on a surface or traces back to the light source
    if (!hit.happened) 
        return;

    // std::cout<<"after eliminating useless photons"<<std::endl;
    
    Vector3f N = normalize(hit.normal);
    Vector3f wo = normalize(-ray.direction);

    // Russian Roulette
    float ksi = get_random_float();
    
    // absorb photon
    if(ksi > RussianRoulette)
        return;

    // reflect photon

    // sample photon bouce direction
    Vector3f reflect_dir = normalize(hit.m->sample(wo, N)); 

    // cast ray from hit to sampled direction
    Ray reflectedPhot = Ray(hit.coords, reflect_dir);
    Intersection reflectedPoint = intersect(reflectedPhot);


    // cannot directly access material type, will result to seg fault.
    // need to find another way to access it
    // if(reflectedPoint.m->m_type == DIFFUSE){
    //     ray = Ray(reflectedPoint.coords, -reflect_dir);
    //     power = power / RussianRoulette;
    //     return;
    // }

    ray = Ray(reflectedPoint.coords, -reflect_dir);
    power = power / RussianRoulette;
    return;
}

void Scene::emitPhotons() const{
    int currNumPhot = 0;
    
    // Object *light;
    // for (uint32_t k = 0; k < objects.size(); ++k)
    //     if (objects[k]->hasEmit()){
    //         light = objects[k];
    //         continue; // since we only have one light source
    //     }
    
    // std::cout<<"inside emitPhotons"<<std::endl;

    Vector3f test = Vector3f(0.0);

    while(currNumPhot < numPhotons){
        Intersection area_p; // to get position on light source
        float pdf_light_area; // no use?
        sampleLight(area_p, pdf_light_area);
        // std::cout<<typeid(area_p.obj).name()<<std::endl;
        // std::cout<<area_p.coords<<std::endl;

        // sample direction
        float x_dir, y_dir, z_dir;
        do{
            x_dir = 2 * get_random_float() - 1;
            y_dir = 2 * get_random_float() - 1;
            z_dir = 2 * get_random_float() - 1;
        } while(x_dir * x_dir + y_dir * y_dir + z_dir * z_dir > 1);
        
        // generate ray from light source and trace photon
        Ray storedPhoton = Ray(area_p.coords, Vector3f(x_dir, y_dir, z_dir));
        Vector3f photonPower = area_p.emit;
        tracePhoton(storedPhoton, photonPower);
        
        // std::cout<<"before storing photons"<<std::endl;

        float dir[3] = {storedPhoton.direction.x, storedPhoton.direction.y, storedPhoton.direction.z};
        float pos[3] = {storedPhoton.origin.x, storedPhoton.origin.y, storedPhoton.origin.z};
        float power[3] = {photonPower.x, photonPower.y, photonPower.z};
        
        // std::cout<<dir<<std::endl;

        // store photons
        causticsMap->store(power, pos, dir);
        globalMap->store(power, pos, dir);

        currNumPhot++;
    }

    std::cout<<"done storing"<<std::endl;

    // scale power
    causticsMap->scale_photon_power(1.0 / numPhotons);
    globalMap->scale_photon_power(1.0 / numPhotons);

    // std::cout<<"no errors up to here"<<std::endl;

    std::cout<<causticsMap->get_num_photons()<<std::endl;

    return;
}

Vector3f Scene::getIrradiance(const Ray &ray) const{
    float irrad[3];
    Intersection its = intersect(ray);

    // if (!its.happened) 
    //     return Vector3f(0);
    // if(its.emit.norm() > 0)
    //     return its.emit;

    float pos[3] = {its.coords.x, its.coords.y, its.coords.z};
    float normal[3] = {-its.normal.x, -its.normal.y, -its.normal.z};
    causticsMap->irradiance_estimate(irrad, pos, normal, 0.1, 100);

    // std::cout<<Vector3f(irrad[0], irrad[1], irrad[2])<<std::endl;


    return Vector3f(irrad[0], irrad[1], irrad[2]);
}


// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray) const {
    Intersection intersection = intersect(ray);
    if (!intersection.happened) {
        //std::cout << ".. ";
        return Vector3f(0.0f);
    }
        

    if (intersection.emit.norm() > 0.f)
         return intersection.m->getEmission();

    const Vector3f pos = intersection.coords;         // position of shading point
    const Vector3f normal = intersection.normal;   // normal at shading point
    const Vector3f obs_dir = -ray.direction;         // observer direction
    Material* m_its = intersection.m;             // material at shading point

    // Direct illumination
    Vector3f L_direct = { 0.0f, 0.0f, 0.0f };

    // sample light sources
    Intersection area_p;
    float pdf_light_area;
    sampleLight(area_p, pdf_light_area);
    //const auto light_sample = scene.sample_light_sources();
    //if (light_sample) {
        const Vector3f light_sample_pos = area_p.coords;               // position of sample point
        Vector3f w_to_light = light_sample_pos - pos;           // shading point to light sample point
        const Vector3f light_sample_dir = w_to_light.normalized();    // light sample point direction
        const Vector3f light_dir = -light_sample_dir;                                   // light direction
        const Vector3f light_sample_normal = area_p.normal;         // normal at sample point

        // check block between shading point and light sample point
        Ray checkerRay = Ray(pos, light_sample_dir);
        const Intersection check_intersection = intersect(checkerRay);
        if (check_intersection.happened && pow((check_intersection.coords - light_sample_pos).norm(),2) < EPSILON) {
            const Vector3f emission = area_p.emit;

            // light source importance sampling
            //const auto light_dir_dot_light_sample_normal = dotProduct(light_dir,light_sample_normal);
            const float pdf_light = (dotProduct(light_dir,light_sample_normal) == 0.0f) ? 
                0.0f : (pow(w_to_light.norm(),2) * pdf_light_area) / fabs(dotProduct(light_dir,light_sample_normal));

            // bsdf importance sampling
            const float pdf_bsdf = m_its->pdf(light_sample_dir, obs_dir, normal);

            // balanced heuristic multiple importance sampling
            const float pdf_sum = pdf_light + pdf_bsdf;
            if (pdf_sum > 0.0f) {
                // std::cout << normal << "??" << light_sample_dir <<"!!";
                L_direct += emission * m_its->eval(light_sample_dir, obs_dir, normal) * fabs(dotProduct(normal, light_sample_dir)) / pdf_sum;
            }
        }
    //}

    // Indirect illumination
    Vector3f L_indirect = { 0.0f, 0.0f, 0.0f };

    // Use Russian Roulette to limit the recursion depth
    if (get_random_float() < RussianRoulette) {
        // sample a direction for indirect illumination
        const Vector3f wi_indir = m_its->sample(obs_dir, normal);

        // bsdf importance sampling
        const float pdf_bsdf = m_its->pdf(wi_indir, obs_dir, normal);
        Ray newRay = Ray(pos, wi_indir);
        if (pdf_bsdf > 0.0f) {
            L_indirect = castRay(newRay) * m_its->eval(wi_indir, obs_dir, normal) * fabs(dotProduct(wi_indir, normal))
                         / (pdf_bsdf * RussianRoulette);
        }
    }

    return  L_direct + L_indirect;

    /* old version castRay */

    // Intersection its = intersect(ray);

    // if (!its.happened) return Vector3f(0.f);
    
    // Vector3f N = normalize(its.normal);
    // Vector3f wo = normalize(-ray.direction);

    // /*
    //     TODO: 
    //     Implement path tracing here.
    // */

    // if(its.emit.norm() > 0)
    //     return its.m->getEmission();

    // // DIRECT ILLUMINATION - contribution from light source
    // // sample light from light source
    // Intersection area_p;
    // float pdf_light_area;
    // sampleLight(area_p, pdf_light_area);

    // // convert to solid angle pdf_omega = ||x' - x||^2/Acos(θ')
    // float light_obj_dist = ((area_p.coords - its.coords).norm()) * ((area_p.coords - its.coords).norm());
    // float area_theta = dotProduct(normalize(its.coords - area_p.coords), normalize(area_p.normal));
    // float pdf_light_solidAngle = light_obj_dist * pdf_light_area / area_theta;
    
    // // cast ray from its to light
    // Vector3f wi_dir = normalize(area_p.coords - its.coords);
    // Ray out_dir = Ray(its.coords, wi_dir);
    // Intersection hitPoint_dir = intersect(out_dir);
    
    // Vector3f L_i = area_p.emit;
    // Vector3f f_r_dir = its.m->eval(wi_dir, wo, N);
    // float theta_dir = dotProduct(wi_dir, N);

    // Vector3f wi = area_p.coords - its.coords;
    // // L_dir = L_i * f_r * cos θ / pdf_light(ω from p to x’)
    // Vector3f L_dir = Vector3f(0.f);
    // if(theta_dir >= cos(M_PI/2) && hitPoint_dir.emit.norm() == area_p.emit.norm())
    // //if(theta_dir >= cos(M_PI/2) && wi.norm() < hitPoint_dir.distance)
    // //&& dotProduct(-wi_dir, normalize(area_p.normal)) > 0
    //     L_dir = L_i * f_r_dir * theta_dir / pdf_light_solidAngle;

    // // INDIRECT ILLUMINATION - contribution from other reflectors
    // Vector3f L_indir = Vector3f(0.f);

    // // Russian Roulette
    // float ksi = get_random_float();
    // if(ksi > RussianRoulette)
    //     return L_dir + L_indir;

    // // Randomly sample the hemisphere toward ωi
    // Vector3f wi_indir = normalize(its.m->sample(wo, N)); 

    // // cast ray from its to sampled direction
    // Ray out_indir = Ray(its.coords, wi_indir);
    // Intersection hitPoint_indir = intersect(out_indir);

    // float pdf_brdf = its.m->pdf(wi_indir, wo, N);
    // Vector3f f_r_indir = its.m->eval(wi_indir, wo, N);
    // float theta_indir = dotProduct(wi_indir, N);
    
    // if(hitPoint_indir.emit.norm() != area_p.emit.norm()) 
    // //if (hitPoint_indir.happened && !hitPoint_indir.obj->hasEmit())
    //     L_indir = castRay(out_indir) * f_r_indir * theta_indir / pdf_brdf / RussianRoulette;

    // return L_indir;
}