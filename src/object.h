#pragma once

#include "ray.h"
#include "vector.h"

struct Object {
    Vector3f albedo;
    Vector3f emittance;
    void setMaterial(Vector3f Albedo, Vector3f Emittance) {
        albedo = Albedo;
        emittance = Emittance;
    }
    virtual Vector3f getNormal(const Vector3f &) const = 0;
    virtual float intersect(const Ray &) const = 0;
};

struct Plane : Object {
    Vector3f position, normal;
    Plane(Vector3f Position, Vector3f Normal, Vector3f baseColor, Vector3f emittance) {
        position = Position;
        normal = Normal;
        this->setMaterial(baseColor, emittance);
    }
    Vector3f getNormal(const Vector3f &x) const { return normal; }
    float intersect(const Ray &ray) const {
        float denominator = Dot(normal, ray.direction);
        if (fabs(denominator) > 1e-6F) {
            float t = Dot(normal, position - ray.origin) / denominator;
            if (t < 1e-6F) return INFINITY;
            return t;
        }
        return INFINITY;
    }
};

struct Sphere : Object {
    float radius, radius2, iradius;
    Vector3f position;
    Sphere(Vector3f Center, float Radius, Vector3f albedo, Vector3f emittance) {
        position = Center;
        radius = Radius;
        radius2 = Radius * Radius;
        this->setMaterial(albedo, emittance);
    }
    Vector3f getNormal(const Vector3f &x) const { return (x - position) / radius; }
    float intersect(const Ray &ray) const {
        Vector3f L = position - ray.origin;
        float tc = Dot(L, ray.direction);
        // if (tc < -1e-6F) return INFINITY;
        float d2 = Dot(L, L) - tc * tc;
        if (d2 > radius2) return INFINITY;
        float t1c = sqrt(radius2 - d2);
        float t = tc - t1c;
        if (t < 1e-6F) {
            t = tc + t1c;
            if (t < 1e-6F) return INFINITY;
            return t;
        }
        return t;
    }
};