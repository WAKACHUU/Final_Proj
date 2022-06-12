#include "Microfacet.hpp"


float Microfacet::distribution(float normal_dot_micro_surface_normal, float roughness_sq) {
    const auto normal_dot_micro_surface_normal_sq = normal_dot_micro_surface_normal * normal_dot_micro_surface_normal;
    auto denominator = normal_dot_micro_surface_normal_sq * (roughness_sq - 1.0f) + 1.0f;
    denominator = M_PI * denominator * denominator;
    return roughness_sq / denominator;
}

Vector3f Microfacet::fresnel_schlick(float micro_surface_normal_dot_ray_out_dir, const Vector3f& f0) {
    return f0 + (Vector3f(1.0f) - f0) * pow(1.0f - micro_surface_normal_dot_ray_out_dir, 5);
}

float Microfacet::geometry(float normal_dot_light_source_dir, float normal_dot_observer_dir, float roughness) {
    return 2.0f / lerp(fabs(2 * normal_dot_light_source_dir * normal_dot_observer_dir), fabs(normal_dot_light_source_dir + normal_dot_observer_dir), roughness);
}

Vector3f Microfacet::sample_micro_surface(const Vector3f& normal, float roughness_sq) {
    const auto r0 = get_random_float();
    const auto r1 = get_random_float();
    const auto theta = std::acos(std::sqrt((1.0f - r0) / ((roughness_sq - 1.0f) * r0 + 1.0f)));
    //float theta = std::acos(std::sqrt(roughness_sq*r0/(1-r0)));
    const auto phi = 2 * M_PI * r1;
	
    const auto local_micro_surface_normal = polar_to_cartesian(theta, phi);
    return toWorld(local_micro_surface_normal, normal);
}

float Microfacet::pdf_micro_surface(float normal_dot_micro_surface_normal, float roughness_sq) {
	// importance sampling on NDF
	const auto normal_dot_micro_surface_normal_abs = fabs(normal_dot_micro_surface_normal);
    return (distribution(normal_dot_micro_surface_normal_abs, roughness_sq) * normal_dot_micro_surface_normal_abs);
}

// Vector3f Microfacet::outward_micro_surface_normal(const Vector3f& ray_source_dir, const Vector3f& ray_out_dir, bool is_same_side, bool is_surface_outward, float ior) {
//     if (is_same_side) {
//         // reflection
//         if (is_surface_outward) {
//             return (ray_source_dir + ray_out_dir).normalized();
//         } else {
//             return -(ray_source_dir + ray_out_dir).normalized();
//         }
//     } else {
//         // refraction
//         if (is_surface_outward) {
//             return -(ray_out_dir + ray_source_dir * ior).normalized();
//         } else {
//             return -(ior * ray_out_dir + ray_source_dir).normalized();
//         }
//     }
// }
Vector3f Microfacet::outward_micro_surface_normal(const Vector3f& wi, const Vector3f& wo, float check_ray, float check_n, float ior) {
    if (check_ray > 0.f) {
        // reflection
        if (check_n > 0.f) {
            //std::cout<<"t1\n";
            //n pointing out
            return (wi + wo).normalized();
        } else {
            //std::cout<<"t2\n";
            //n pointing in
            return -(wi + wo).normalized();
        }
    } else {
        // refraction
        if (check_n > 0.f) {
            //std::cout<<"t3\n";
            return -(wo + wi * ior).normalized();
        } else {
            //std::cout<<"t4\n";
            return -(ior * wo + wi).normalized();
        }
    }
}

float Microfacet::reflect_jacobian(float micro_surface_normal_dot_ray_out_dir) {
    return micro_surface_normal_dot_ray_out_dir == 0.0f ? 0.0f : 1.0f / (4 * fabs(micro_surface_normal_dot_ray_out_dir));
}

float Microfacet::refract_jacobian(float micro_surface_normal_dot_ray_source_dir, float ior_in, float micro_surface_normal_dot_ray_out_dir, float ior_out) {
    auto denominator = ior_in * micro_surface_normal_dot_ray_source_dir + ior_out * micro_surface_normal_dot_ray_out_dir;
    denominator *= denominator;
    return denominator == 0.0f ? 0.0f : (ior_out * ior_out * fabs(micro_surface_normal_dot_ray_out_dir)) / denominator;
}