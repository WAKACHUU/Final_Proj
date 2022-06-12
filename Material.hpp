//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"
#include "Microfacet.hpp"

using std::swap;

enum MaterialType {DIFFUSE, MIRROR, TR};

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

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
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

    inline float sampleGlass(const Vector3f &wi, const Vector3f &wo, const Vector3f &N, float* pdf);

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
        {
            // randomly choose a micro surface
            const auto micro_surface_normal = Microfacet::sample_micro_surface(N, roughness*roughness);
            const auto obs_dir  = -wi;

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
            if(dotProduct(wo, N) > 0.0f)
                return 1.0f;
            else
                return 0.0f;
            break;
        }
        
        case TR:
        {
            // mac
            const auto cosine_i = dotProduct(wi, N);
            const auto cosine_o = dotProduct(wo, N);

            const auto check_ray_dir = cosine_i * cosine_o;

            // mic
            const auto micro_surface_normal = Microfacet::outward_micro_surface_normal(wi, wo, check_ray_dir, cosine_o, ior);
            const auto cosine_i_mf = dotProduct(wi, micro_surface_normal);
            const auto cosine_o_mf = dotProduct(wo, micro_surface_normal);

            // extreme case
            if (check_ray_dir == 0.0f)
                return 0.0f;

            const auto obs_dir = -wo;

            float f;
            fresnel(obs_dir, micro_surface_normal, ior, f);

            const auto pdf_micro_surface = Microfacet::pdf_micro_surface(dotProduct(micro_surface_normal, N), roughness*roughness);



            if (check_ray_dir > 0.f) {
                // reflection
                const auto jacobian = Microfacet::reflect_jacobian(cosine_o_mf);
                return pdf_micro_surface * f * jacobian;
            } else {
                // refraction
                const auto ior_in = cosine_i < 0.0f ? ior : 1.0f;
                const auto ior_out = cosine_o < 0.0f ? ior : 1.0f;
                const auto jacobian = Microfacet::refract_jacobian(cosine_i_mf, ior_in, cosine_o_mf, ior_out);
                return pdf_micro_surface * (1.0f - f) * jacobian;
            }
        }
        
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
            if(cosalpha > 0.0f){
                float kr;
                fresnel(wi, N, ior, kr);
                return kr * Vector3f(1.0f / cosalpha);
            } else
                return Vector3f(0.0f);

            break;
        }

        case TR: 
        // {
        //     float cosalpha = dotProduct(N, wo);
        //     if(cosalpha < 0.001f){
        //         std::cout << ".. ";
        //         float kr;
        //         fresnel(wi, N, 1.5f, kr);
        //         return (1.0f-kr) * Vector3f(1.0f / -cosalpha);
        //     }
        // }
        {
            // mac
            const auto cosine_i = dotProduct(wi, N);
            const auto cosine_o = dotProduct(wo, N);
            const auto check_ray_dir = cosine_i * cosine_o;

            // mic
            const auto micro_surface_normal = Microfacet::outward_micro_surface_normal(wi, wo, check_ray_dir, cosine_o, ior);
            const auto cosine_i_mf = dotProduct(wi, micro_surface_normal);
            const auto cosine_o_mf = dotProduct(wo, micro_surface_normal);


            if (check_ray_dir == 0.0f)
                return 0.0f;

            
            // const auto is_same_side = check_ray_dir > 0.0f;
            // const auto is_surface_outward = normal_dot_ray_out_dir > 0.0f;



            //const auto normal_dot_micro_surface_normal = dotProduct(micro_surface_normal, N);


            const auto D = Microfacet::distribution(dotProduct(micro_surface_normal, N), roughness*roughness);
            const auto G = Microfacet::geometry(cosine_i_mf, cosine_o_mf, roughness);

            const auto obs_dir = -wo;
            float F;
            fresnel(obs_dir, micro_surface_normal, ior, F);

            if (check_ray_dir > 0.f) {
                // reflection
                //std::cout << "eval refl" << D * F * G / 4.0f << std::endl;
                return D * F * G / 4.0f;   // Cook–Torrance Specular (original denominator is merged into G for Smith-Joint approximation)
            } else {
                // refraction
                const auto ior_in = cosine_i < 0.0f ? ior : 1.0f;
                const auto ior_out = cosine_o < 0.0f ? ior : 1.0f;

                // std::cout << "eval: refr:" << Microfacet::refract_jacobian(micro_surface_normal_dot_ray_source_dir, ior_in, micro_surface_normal_dot_ray_out_dir, ior_out)
                //     << " " << fabs(micro_surface_normal_dot_ray_source_dir) <<" "<< D * (1.0f - F) * G <<std::endl;

                // Transmission (original denominator is merged into G for Smith-Joint approximation)
                return Microfacet::refract_jacobian(cosine_i_mf, ior_in, cosine_o_mf, ior_out)
                    * fabs(cosine_i_mf) * D * (1.0f - F) * G;
            }
        }
    }

}
float Material::sampleGlass(const Vector3f &wi, const Vector3f &wo, const Vector3f &N, float* pdf){
    
    // float R0 = (ior-1.0)*(ior-1.0)/((ior + 1.0)*(ior + 1.0));
    // float cosTheta = clamp(-1, 1, dotProduct(wi, N));
    // float f = 1 - fabs(cosTheta);
    // float g = ((f * f) * (f * f)) * f;
    // float fresnel_coe = R0 + (1.0 - R0)*g;

    // bool entering = cosTheta > 0;
    // float ei = 1.f, et = ior;
    // if (!entering) {
    //     swap(ei, et);
    //     cosTheta = -cosTheta;  // be careful here, want cosTheta to be
    //                             // positive for everything below
    // }

    // float inveta = et / ei;
    // float inveta2 = inveta * inveta;

    // if (refract(wi, N, ior) == 0) {
    //     // total internal reflection; always reflect
    //     *pdf = 1.0;
    //     // reflect(wi, N);
    //     return Vector3f(1.0f / cosTheta);
    // }

    // if ((double)(std::rand()) / RAND_MAX < fresnel_coe) {
    //     *pdf = fresnel_coe;
    //     // reflect(wo, wi);
    //     return (fresnel_coe / cosTheta) * Vector3f(1.f);
    // } else {
    //     // refraction ray has already been computed
    //     float one_minus_fresnel = 1.0f - fresnel_coe;
    //     *pdf = one_minus_fresnel;
    //     return (one_minus_fresnel * inveta2 / cosTheta) * Vector3f(1.f);
    // }

    // refract(wi, N, ior);
    
    Vector3f sampledDir = refract(wi, N, ior);
    if(sampledDir == 0) return 0;
    float cosTheta = dotProduct(sampledDir, -N);

    float ft;
    fresnel(wi, N, ior, ft);
    ft = 1 - ft;

    // (etaI * etaI) / (etaT * etaT);
    return (1 / ior * ior) * ft / cosTheta;

}

#endif //RAYTRACING_MATERIAL_H
