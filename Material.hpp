//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"
#include "Microfacet.hpp"

enum MaterialType { DIFFUSE, MIRROR, TR};

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }


public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    float roughness;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wi, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            
            break;
        }

        case MIRROR:
        {
            Vector3f localRay = reflect(wi, N);
            return toWorld(localRay, N);
            break;
        }

        case TR:
        {
            // randomly choose a micro surface
            const Vector3f micro_surface_normal = Microfacet::sample_micro_surface(N, roughness*roughness);
            const Vector3f obs_dir  = -wi;

            float f;
            fresnel(obs_dir, micro_surface_normal, ior, f);

            // trace back
            if (get_random_float() < f) {
                // reflection
                return reflect(obs_dir, micro_surface_normal);
            } else {
                // refraction(transmission)
                return refract(obs_dir, micro_surface_normal, ior);
            }
        }
        // {
        //     Vector3f localRay = refract(wi, N, 1.5f);

        //     Vector3f test = Vector3f(1,1,1).normalized();
        //     Vector3f testn = Vector3f(0,1,0);
        //     std::cout << "??";
        //     std::cout << refract(test, testn, 1.5f).x << refract(test, testn, 1.5f).y << refract(test, testn, 1.5f).z;
        //     //std::cout << dotProduct(localRay,N);
        //     return localRay;

        //     break;
        // }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case MIRROR:
        {
            if(dotProduct(wo, N) > 0.001f)
                return 1.0f;
            else
                return 0.0f;
            break;
        }
        case TR:
        {
            const float check_ray_dir = dotProduct(wi, N) * dotProduct(wo, N);

            
            if (check_ray_dir == 0.0f) {
                return 0.0f;
            }
                
            const Vector3f obs_dir = -wo;

            const bool is_same_side = check_ray_dir > 0.0f;
            const bool is_surface_outward = dotProduct(wo, N) > 0.0f;

            const Vector3f micro_surface_normal = Microfacet::outward_micro_surface_normal(wi, wo, is_same_side, is_surface_outward, ior);
            float f;
            fresnel(obs_dir, micro_surface_normal, ior, f);

            const float normal_dot_micro_surface_normal = dotProduct(micro_surface_normal, N);
            const float pdf_micro_surface = Microfacet::pdf_micro_surface(normal_dot_micro_surface_normal, roughness*roughness);

            if (is_same_side) {
                // reflection
                const auto jacobian = Microfacet::reflect_jacobian(dotProduct(wo, micro_surface_normal));
                return pdf_micro_surface * f * jacobian;
            } else {
                // refraction
                //determine in and out
                const auto ior_in = dotProduct(wi, N) < 0.0f ? ior : 1.0f;
                const auto ior_out = dotProduct(wo, N) < 0.0f ? ior : 1.0f;

                const auto jacobian = Microfacet::refract_jacobian(dotProduct(wi, micro_surface_normal), ior_in, dotProduct(wo, micro_surface_normal), ior_out);
                return pdf_micro_surface * (1.0f - f) * jacobian;
            }
        }
        // {
        //     if (dotProduct(wo, N) < 0.001f) {
        //         return 1.0f;
        //     }
        //     else {
        //         return 0.0f;
        //     }
        // }
    }
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MIRROR:
        {
            float cosalpha = dotProduct(N, wo);
            if(cosalpha > 0.001f){
                float kr;
                fresnel(wi, N, ior, kr);
                return kr * Vector3f(1.0f / cosalpha);
            }
            else {
                return Vector3f(0.0f);
            }
        }
        case TR: 
        {
            const float normal_dot_wi = dotProduct(wi, N);
            const float normal_dot_wo = dotProduct(wo, N);
            const float check_ray_dir = normal_dot_wi * normal_dot_wo;
            if (check_ray_dir == 0.0f)
                return 0.0f;

            const Vector3f obs_dir = -wo;
            const bool is_same_side = check_ray_dir > 0.0f;
            const bool is_surface_outward = normal_dot_wo > 0.0f;

            const Vector3f micro_surface_normal = Microfacet::outward_micro_surface_normal(wi, wo, is_same_side, is_surface_outward, ior);


            const auto D = Microfacet::distribution(dotProduct(micro_surface_normal, N), roughness*roughness);
            const auto G = Microfacet::geometry(dotProduct(wi, micro_surface_normal), 
                                                dotProduct(wo, micro_surface_normal), 
                                                roughness);
            float F;
            fresnel(obs_dir, micro_surface_normal, ior, F);

            if (is_same_side) {
                // reflection
                //std::cout << "eval refl" << D * F * G / 4.0f << std::endl;
                return D * F * G / 4.0f;   // Cook–Torrance Specular (original denominator is merged into G for Smith-Joint approximation)
            } else {
                // refraction
                const auto ior_in = normal_dot_wi < 0.0f ? ior : 1.0f;
                const auto ior_out = normal_dot_wo < 0.0f ? ior : 1.0f;

                // std::cout << "eval: refr:" << Microfacet::refract_jacobian(dotProduct(wi, micro_surface_normal), ior_in, dotProduct(wo, micro_surface_normal), ior_out)
                //     << " " << fabs(dotProduct(wi, micro_surface_normal)) <<" "<< D * (1.0f - F) * G <<std::endl;

                // Transmission (original denominator is merged into G for Smith-Joint approximation)
                return Microfacet::refract_jacobian(dotProduct(wi, micro_surface_normal), ior_in, dotProduct(wo, micro_surface_normal), ior_out)
                    * fabs(dotProduct(wi, micro_surface_normal)) * D * (1.0f - F) * G;
            }
        }
        // {
        //     float cosalpha = dotProduct(N, wo);
        //     if(cosalpha < 0.001f){
        //         std::cout << ".. ";
        //         float kr;
        //         fresnel(wi, N, 1.5f, kr);
        //         return (1.0f-kr) * Vector3f(1.0f / -cosalpha);
        //     }
        // }
    }
}

#endif //RAYTRACING_MATERIAL_H
