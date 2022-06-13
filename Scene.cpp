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

void Scene::tracePhoton(Ray &ray, Vector3f &power, bool hitSpecular, int depth) const{
    
    // std::cout<<"inside trace photons"<<std::endl;

    Intersection hit = intersect(ray);

    // photon doesn't land on a surface or traces back to the light source
    if (!hit.happened) 
        return;

    // std::cout<<"after eliminating useless photons"<<std::endl;

    // photon ray hitting a diffuse surface
    float power_arr[3] = {power.x, power.y, power.z};
    float pos[3] = {hit.coords.x, hit.coords.y, hit.coords.z};
    float dir[3] = {ray.direction.x, ray.direction.y, ray.direction.z};
    if(hit.m->getType() == DIFFUSE){
        if(hitSpecular){ // hit specular/glossy at least once
            causticsMap->store(power_arr, pos, dir);
            return;
        } 

        globalMap->store(power_arr, pos, dir); 

        // Russian Roulette to decide if continue to bounce within the scene or absorb the photon
        // absorb photon = stop the photon tracing for the current sample
        float ksi = get_random_float();
        if(ksi > RussianRoulette)
            return;
    } 
    
    // continue to bounce: applies to all material
    if(depth > 20) // there will be a segfault if don't limit the number bounces of the photon ray
        return;

    Vector3f N = normalize(hit.normal);
    Vector3f w_in = normalize(-ray.direction);
    Vector3f w_new = normalize(hit.m->sample(w_in, N)); 
    Vector3f f = hit.m->eval(w_new, w_in, N); 
    float pdf = hit.m->pdf(w_in, w_new, N);
    float cosTheta = dotProduct(w_new, N);

    // Ray out_indir = Ray(its.coords, wi_indir);
    // Intersection hitPoint_indir = intersect(out_indir);
    
    Vector3f R = (f * cosTheta) / pdf;
    if(hit.m->getType() == TR) R = R / 2;

    Ray newBounce = Ray(hit.coords, w_new);
    // TODO: figure out how to scale the power, the 10 here is just a temp value
    Vector3f newPower = power * 0.1 / RussianRoulette;
    // Vector3f origPower = power;
    if(hit.m->getType() == MIRROR || hit.m->getType() == TR) 
        tracePhoton(newBounce, power, true, ++depth); // don't scale power specular?
    else tracePhoton(newBounce, newPower, false, ++depth);
    
    return;
}

void Scene::photonMapping() const{    
    
    // std::cout<<"inside emitPhotons"<<std::endl;

    Vector3f test = Vector3f(0.0);

    int currNumPhotons = 0;
    while(currNumPhotons < numPhotons){
        Intersection area_p; // to get position on light source
        float pdf_light_area; // no use?
        sampleLight(area_p, pdf_light_area);

        // std::cout<<area_p.normal<<std::endl; // (0, -1, 0)

        // sample direction from sampled point on light source
        float u = get_random_float();
        float v = 2 * M_PI * get_random_float();
        Vector3f sampledDir = Vector3f(cos(v) * sqrt(u), sin(v) * sqrt(u), sqrt(1 - u));

        // calculate power based on sampled direction
        float cosTheta = dotProduct(normalize(sampledDir), area_p.normal);
        float pdf_light_dir = cosTheta / M_PI;
        Vector3f power_vec = (area_p.emit * cosTheta) / (pdf_light_area * pdf_light_dir);

        // generate ray from light source and start tracing a photon
        Ray storedPhoton = Ray(area_p.coords,  sampledDir); // Vector3f(x_dir, y_dir, z_dir)
        tracePhoton(storedPhoton, power_vec, false, 0);
        
        currNumPhotons++;
    }

    std::cout<<"done storing"<<std::endl;

    causticsMap->balance();
    globalMap->balance();
    causticsMap->scale_photon_power(1.0 / causticsMap->get_num_photons());
    globalMap->scale_photon_power(1.0 / numPhotons);

    // std::cout<<"no errors up to here"<<std::endl;
    std::cout<<causticsMap->get_num_photons()<<std::endl;
    std::cout<<globalMap->get_num_photons()<<std::endl;

    return;
}


// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray) const {

    Intersection intersection = intersect(ray);

    Vector3f N = normalize(intersection.normal);
    Vector3f wo = normalize(-ray.direction);

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
// //
//     //const auto light_sample = scene.sample_light_sources();
//     //if (light_sample) {
//         const Vector3f light_sample_pos = area_p.coords;               // position of sample point
//         Vector3f w_to_light = light_sample_pos - pos;           // shading point to light sample point
//         const Vector3f light_sample_dir = w_to_light.normalized();    // light sample point direction
//         const Vector3f light_dir = -light_sample_dir;                                   // light direction
//         const Vector3f light_sample_normal = area_p.normal;         // normal at sample point

//         // check block between shading point and light sample point
//         Ray checkerRay = Ray(pos, light_sample_dir);
//         const Intersection check_intersection = intersect(checkerRay);
//         if (check_intersection.happened && pow((check_intersection.coords - light_sample_pos).norm(),2) < EPSILON) {
//             const Vector3f emission = area_p.emit;

//             // light source importance sampling
//             //const auto light_dir_dot_light_sample_normal = dotProduct(light_dir,light_sample_normal);
//             const float pdf_light = (dotProduct(light_dir,light_sample_normal) == 0.0f) ? 
//                 0.0f : (pow(w_to_light.norm(),2) * pdf_light_area) / fabs(dotProduct(light_dir,light_sample_normal));

//             // bsdf importance sampling
//             const float pdf_bsdf = m_its->pdf(light_sample_dir, obs_dir, normal);

//             // balanced heuristic multiple importance sampling
//             const float pdf_sum = pdf_light + pdf_bsdf;
//             if (pdf_sum > 0.0f) {
//                 // std::cout << normal << "??" << light_sample_dir <<"!!";
//                 L_direct += emission * m_its->eval(light_sample_dir, obs_dir, normal) * fabs(dotProduct(normal, light_sample_dir)) / pdf_sum;
//             }
//         }
//     //}

//     // Indirect illumination
//     Vector3f L_indirect = { 0.0f, 0.0f, 0.0f };

//     // Use Russian Roulette to limit the recursion depth
//     if (get_random_float() < RussianRoulette) {
//         // sample a direction for indirect illumination
//         const Vector3f wi_indir = m_its->sample(obs_dir, normal);

//         // bsdf importance sampling
//         const float pdf_bsdf = m_its->pdf(wi_indir, obs_dir, normal);
//         Ray newRay = Ray(pos, wi_indir);
//         if (pdf_bsdf > 0.0f) {
//             L_indirect = castRay(newRay) * m_its->eval(wi_indir, obs_dir, normal) * fabs(dotProduct(wi_indir, normal))
//                          / (pdf_bsdf * RussianRoulette);
//         }
//     }
//
    // if(depth > 1)
    //     return  L_direct + L_indirect;

    // INDIRECT ILLUMINATION - GLOBAL PHOTON MAP
    Vector3f L_gm = Vector3f(0.0), L_c = Vector3f(0.0);
    if(intersection.m->getType() == DIFFUSE){
        // std::cout<<"inside indirect"<<std::endl;

        // sample hemisphere toward ωi
        Vector3f wi_pm = normalize(intersection.m->sample(wo, N)); 
        Vector3f f_c = intersection.m->eval(wi_pm, wo, N); // not sure about brdf
        Vector3f f_gm = f_c; // probs have to change this later
        Ray out_gm = Ray(intersection.coords, wi_pm);
        Intersection hitPoint_gm = intersect(out_gm);
        // std::cout<<hitPoint_gm.coords<<std::endl;

        // TODO: don't know how to get material for hitPoint_gm since it's always giving segfault
        // need material to calculate brdf at the new bounce?
        
        if(hitPoint_gm.emit.norm() != area_p.emit.norm()){

            // global PM
            // bounce twice
            float pos_gm[3] = {hitPoint_gm.coords.x, hitPoint_gm.coords.y, hitPoint_gm.coords.z};
            float normal_gm[3] = {hitPoint_gm.normal.x, hitPoint_gm.normal.y, hitPoint_gm.normal.z};
            float irrad_gm[3] = {area_p.emit.x, area_p.emit.y, area_p.emit.z};
            globalMap->irradiance_estimate(irrad_gm, pos_gm, normal_gm, 1.5, 100);
            Vector3f irrad_gm_vec = Vector3f(irrad_gm[0], irrad_gm[1], irrad_gm[2]);
            // std::cout<<irrad_gm_vec<<std::endl;
            L_gm = f_gm * irrad_gm_vec; // result is currently too bright if no scaling involved in photon tracing

            // caustics PM
            float pos_c[3] = {intersection.coords.x, intersection.coords.y, intersection.coords.z};
            float normal_c[3] = {intersection.normal.x, intersection.normal.y, intersection.normal.z};
            float irrad_c[3] = {area_p.emit.x, area_p.emit.y, area_p.emit.z};
            causticsMap->irradiance_estimate(irrad_c, pos_c, normal_c, 1.5, 150);
            Vector3f irrad_c_vec = Vector3f(irrad_c[0], irrad_c[1], irrad_c[2]);
            // std::cout<<irrad_c_vec<<std::endl;
            L_c = f_c * irrad_c_vec; // result is currently too bright if no scaling involved in photon tracing
        }
    }

    // return L_gm + L_dir + L_c;
    //return L_direct + L_indirect + L_gm + L_c;
    // return L_dir + L_indir;
    return L_c;
    // return L_gm;
    // return L_indir;
    // return L_dir;

    // Intersection its = intersect(ray);

    // if (!its.happened) return Vector3f(0.f);
    
    // Vector3f N = normalize(its.normal);
    // Vector3f wo = normalize(-ray.direction);

    // // Le - light emission
    // if(its.emit.norm() > 0)
    //     return its.m->getEmission();

    // // DIRECT ILLUMINATION - contribution from light source to diffuse surfaces
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

    // float f_dir_trans;
    // if(its.m->getType() == REFRACTIVE){
    //     f_dir_trans = its.m->sampleGlass(wo, wi, N, &pdf_light_solidAngle);
    // }
    
    // // float pdf = 0.;
    // // if(its.m->getType() == GLASS){
    // //         its.m->randGlass = get_random_float();
    // //         f_r_dir = its.m->sampleGlass(wi, wo, N, &pdf);
    // //     }

    // // L_dir = L_i * f_r * cos θ / pdf_light(ω from p to x’)
    // Vector3f L_dir = Vector3f(0.f);
    // if(theta_dir >= cos(M_PI/2) && hitPoint_dir.emit.norm() == area_p.emit.norm()){
    //     // if(its.m->getType() == GLASS){
    //     //    float *pdf;
    //     //    f_r_dir = its.m->sampleGlass(wi_dir, wo, N, pdf);
    //     // }

    //     if(its.m->getType() == DIFFUSE)
    //         L_dir = L_i * f_r_dir * theta_dir / pdf_light_solidAngle;

    //     if(its.m->getType() == REFRACTIVE)
    //         L_dir = L_i * f_dir_trans;

    //     // if(its.m->getType() == GLASS)
    //     //     if(pdf == pdf_light_area)
    //     //         L_dir = L_i * f_r_dir * theta_dir / pdf_light_solidAngle;
            
    //     // TODO: continue to bounce if hit specular/glossy surface
    // }

    // return L_dir;

    // // INDIRECT ILLUMINATION - contribution from other reflectors
    // Vector3f L_indir = Vector3f(0.f);

    // // Russian Roulette
    // float ksi = get_random_float();
    // if(ksi < RussianRoulette){
    //     // Randomly sample the hemisphere toward ωi
    //     Vector3f wi_indir = normalize(its.m->sample(wo, N)); 

    //     // cast ray from its to sampled direction
    //     Ray out_indir = Ray(its.coords, wi_indir);
    //     Intersection hitPoint_indir = intersect(out_indir);

    //     float pdf_brdf = its.m->pdf(wi_indir, wo, N);
    //     Vector3f f_r_indir = its.m->eval(wi_indir, wo, N);
    //     float theta_indir = dotProduct(wi_indir, N);
        
    //     if(hitPoint_indir.emit.norm() != area_p.emit.norm()) {
    //         // if(its.m->getType() == GLASS)
    //         //     L_indir = castRay(out_indir) * f_r_indir / pdf_brdf / RussianRoulette;
    //         // // if(its.m->getType() == DIFFUSE)
    //         // else 
    //         if(its.m->getType() == GLASS){
    //             its.m->randGlass = get_random_float();
                
    //         }

    //         L_indir = castRay(out_indir) * f_r_indir * theta_indir / pdf_brdf / RussianRoulette;
            
    //     }
    // }

    // // INDIRECT ILLUMINATION - GLOBAL PHOTON MAP
    // Vector3f L_gm = Vector3f(0.0), L_c = Vector3f(0.0);
    // if(its.m->getType() == DIFFUSE){
    //     // std::cout<<"inside indirect"<<std::endl;

    //     // sample hemisphere toward ωi
    //     Vector3f wi_pm = normalize(its.m->sample(wo, N)); 
    //     Vector3f f_c = its.m->eval(wi_pm, wo, N); // not sure about brdf
    //     Vector3f f_gm = f_c; // probs have to change this later
    //     Ray out_gm = Ray(its.coords, wi_pm);
    //     Intersection hitPoint_gm = intersect(out_gm);
    //     // std::cout<<hitPoint_gm.coords<<std::endl;

    //     // TODO: don't know how to get material for hitPoint_gm since it's always giving segfault
    //     // need material to calculate brdf at the new bounce?
        
    //     if(hitPoint_gm.emit.norm() != area_p.emit.norm()){

    //         // global PM

    //         float pos_gm[3] = {hitPoint_gm.coords.x, hitPoint_gm.coords.y, hitPoint_gm.coords.z};
    //         float normal_gm[3] = {hitPoint_gm.normal.x, hitPoint_gm.normal.y, hitPoint_gm.normal.z};
    //         float irrad_gm[3] = {area_p.emit.x, area_p.emit.y, area_p.emit.z};
    //         globalMap->irradiance_estimate(irrad_gm, pos_gm, normal_gm, 20, 500);
    //         Vector3f irrad_gm_vec = Vector3f(irrad_gm[0], irrad_gm[1], irrad_gm[2]);
    //         // std::cout<<irrad_gm_vec<<std::endl;
    //         L_gm = f_gm * irrad_gm_vec; // result is currently too bright if no scaling involved in photon tracing

    //         // caustics PM
    //         float pos_c[3] = {its.coords.x, its.coords.y, its.coords.z};
    //         float normal_c[3] = {its.normal.x, its.normal.y, its.normal.z};
    //         float irrad_c[3] = {area_p.emit.x, area_p.emit.y, area_p.emit.z};
    //         causticsMap->irradiance_estimate(irrad_c, pos_c, normal_c, 50, 100);
    //         Vector3f irrad_c_vec = Vector3f(irrad_c[0], irrad_c[1], irrad_c[2]);
    //         // std::cout<<irrad_c_vec<<std::endl;
    //         L_c = f_c * irrad_c_vec; // result is currently too bright if no scaling involved in photon tracing
    //     }
    // }

    // // return L_gm + L_dir + L_c;
    // // return L_dir + L_indir + L_gm + L_c;
    // // return L_dir + L_indir;
    // // return L_c;
    // // return L_gm;
    // // return L_indir;
    // // return L_dir;

    // return L_dir + L_indir;
    
}