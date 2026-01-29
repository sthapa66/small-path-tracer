#pragma once

#include "vector.h"

struct Ray {
    Vector3f origin, direction;

    Ray() : origin(Vector3f()), direction(Vector3f()) {}

    Ray(Vector3f origin, Vector3f direction) : origin(origin), direction(direction) {}

    Vector3f operator()(const float &t) const {
        return origin + direction * t;
    }
};