#include "Microfacet.hpp"


float Microfacet::distribution(float normal_dot_micro_surface_normal, float roughness_sq) {
    const auto normal_dot_micro_surface_normal_sq = normal_dot_micro_surface_normal * normal_dot_micro_surface_normal;
    auto denominator = normal_dot_micro_surface_normal_sq * (roughness_sq - 1.0f) + 1.0f;
    denominator = M_PI * denominator * denominator;
    return roughness_sq / denominator;
}


float Microfacet::geometry(float normal_dot_light_source_dir, float normal_dot_observer_dir, float roughness) {
    return 2.0f / lerp(fabs(2 * normal_dot_light_source_dir * normal_dot_observer_dir), fabs(normal_dot_light_source_dir + normal_dot_observer_dir), roughness);
}

Vector3f Microfacet::sample_micro_surface(const Vector3f& normal, float roughness_sq) {
    const auto r0 = get_random_float();
    const auto r1 = get_random_float();
    const auto theta = std::acos(std::sqrt((1.0f - r0) / ((roughness_sq - 1.0f) * r0 + 1.0f)));
    const auto phi = 2 * M_PI * r1;
	
    const auto local_micro_surface_normal = polar_transform(theta, phi);
    return toWorld(local_micro_surface_normal, normal);
}

float Microfacet::pdf_micro_surface(float normal_dot_micro_surface_normal, float roughness_sq) {
	// importance sampling on NDF
	const auto normal_dot_micro_surface_normal_abs = fabs(normal_dot_micro_surface_normal);
    return (distribution(normal_dot_micro_surface_normal_abs, roughness_sq) * normal_dot_micro_surface_normal_abs);
}

Vector3f Microfacet::outward_micro_surface_normal(const Vector3f& wi, const Vector3f& wo, bool is_same_side, bool is_surface_outward, float ior) {
    if (is_same_side) {
        // reflection
        if (is_surface_outward) {
            return (wi + wo).normalized();
        } else {
            return -(wi + wo).normalized();
        }
    } else {
        // refraction
        if (is_surface_outward) {
            return -(wo + wi * ior).normalized();
        } else {
            return -(ior * wo + wi).normalized();
        }
    }
}

float Microfacet::reflect_jacobian(float micro_surface_normal_dot_wo) {
    return micro_surface_normal_dot_wo == 0.0f ? 0.0f : 1.0f / (4 * fabs(micro_surface_normal_dot_wo));
}

float Microfacet::refract_jacobian(float micro_surface_normal_dot_wi, float ior_in, float micro_surface_normal_dot_wo, float ior_out) {
    auto denominator = ior_in * micro_surface_normal_dot_wi + ior_out * micro_surface_normal_dot_wo;
    denominator *= denominator;
    return denominator == 0.0f ? 0.0f : (ior_out * ior_out * fabs(micro_surface_normal_dot_wo)) / denominator;
}
/////////////////////////////