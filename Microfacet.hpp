#ifndef RAYTRACING_MICROFACET_H
#define RAYTRACING_MICROFACET_H

#include "Vector.hpp"
#include "global.hpp"

////////////////////////


////-------------


inline float lerp(float x, float y, float t) { return x * (1.0f - t) + y * t; }

inline Vector3f toWorld(const Vector3f &a, const Vector3f &N){
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

// transform polar coordinates to cartesian coordinates
inline Vector3f polar_to_cartesian(float theta, float phi) {
	return {
        std::sin(theta) * std::cos(phi),
        std::sin(theta) * std::sin(phi),
         std::cos(theta)
	};
}

class Microfacet {
public:
    // Generalized-Trowbridge-Reitz Normal Distribution Function (GTR-NDF) when γ = 2.
    static float distribution(float normal_dot_micro_surface_normal, float roughness_sq);
    // Fresnel-Schlick approximation
    static Vector3f fresnel_schlick(float micro_surface_normal_dot_ray_out_dir, const Vector3f& f0);
    // Smith-Joint Approximation from Respawn Entertainment.
    static float geometry(float normal_dot_light_source_dir, float normal_dot_observer_dir, float roughness_sq);

	// Sample a micro-surface under the distribution function and calculate its surface normal.
    static Vector3f sample_micro_surface(const Vector3f& normal, float roughness_sq);
    // Probability distribution function for importance sampling on GTR-NDF
    static float pdf_micro_surface(float normal_dot_micro_surface_normal, float roughness_sq);
    // Calculate the outward micro-surface normal vector
    static Vector3f outward_micro_surface_normal(const Vector3f& ray_source_dir, const Vector3f& ray_out_dir,
                                              float check_ray, float check_n, float ior);

    // The absolute value of the determinant of the Jacobian matrix for the transformation
    // between micro-surface normal and reflected ray.
    static float reflect_jacobian(float micro_surface_normal_dot_ray_out_dir);
    // The absolute value of the determinant of the Jacobian matrix for the transformation
    // between micro-surface normal and refracted ray.
    static float refract_jacobian(float micro_surface_normal_dot_ray_source_dir, float ior_in,
                                  float micro_surface_normal_dot_ray_out_dir, float ior_out);
};


#endif