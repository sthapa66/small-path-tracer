#pragma once

#include "vector.h"

struct Material {
    Vector3f albedo;
    float roughness;
    float metallic;
    Vector3f emittance;
    Material() : albedo(Vector3f(0.8)), roughness(1), metallic(0), emittance() {}
    Material(Vector3f Albedo, float Roughness, float Metallic, Vector3f Emittance) {
        albedo = Albedo;
        roughness = Roughness;
        metallic = Metallic;
        emittance = Emittance;
    }
};