#include <math.h>

#include <ctime>
#include <iostream>
#include <random>
#include <vector>

#include "atmosphere.h"
#include "object.h"
#include "ray.h"
#include "vector.h"

std::mt19937 mersenneTwister;
std::uniform_real_distribution<float> uniform;
#define random uniform(mersenneTwister)

template <typename T>
inline T Saturate(T x) {
    return std::min((T)1, std::max((T)0, x));
}

struct Intersection {
    float t;
    Object *object;
    Intersection() : object(nullptr), t(INFINITY) {}
    Intersection(Object *object, float t) : object(object), t(t) {}
};

struct Scene {
    std::vector<Object *> objects;
    void add(Object *object) { objects.push_back(object); }
    Intersection intersect(const Ray &ray) const {
        Intersection intersection;
        for (auto iter = objects.begin(); iter != objects.end(); iter++) {
            float t = (*iter)->intersect(ray);
            if (t < intersection.t) {
                intersection.object = *iter;
                intersection.t = t;
            }
        }
        return intersection;
    }
};

Vector3f UniformHemisphereSample(Vector3f n) {
    float r = sqrt(random) * 1.57079632679F;
    float t = random * 6.28318530718F;
    float s = sin(r);
    float x = s * cos(t), y = s * sin(t), z = cos(r);
    Vector3f T = Normalize(Cross(n, Vector3f(0, 1, 1)));
    Vector3f B = Cross(n, T);
    return Normalize(T * x + B * y + n * z);
}


Vector3f PathTrace(const Ray &ray, const Scene &scene, const int &depth) {
    if (depth > 4) return Vector3f(0);

    Intersection intersection = scene.intersect(ray);

    if (intersection.object == nullptr) {
        float gr = 6360e3;
        float ar = 6420e3;
        float HR = 7994;
        float HM = 1200;
        Vector3f betaR(0.00000530482236632, 0.0000123952627965, 0.0000302618720618);
        Vector3f betaM(0.000021);
        Vector3f sd = Normalize(Vector3f(0, 1, -2));
        Atmosphere atmosphere(gr, ar, HR, HM, betaR, betaM, sd);
        Ray aViewRay;
        aViewRay.origin = Vector3f(0, gr + 500.0f, 0);
        aViewRay.direction = ray.direction;
        Vector3f sky = ComputeAtmosphereicScattering(aViewRay, atmosphere);
        return Max(sky, 0.0F);
    }
    Vector3f emittance = intersection.object->emittance;
    Vector3f hitPoint = ray.origin + ray.direction * intersection.t;
    Vector3f normal = intersection.object->getNormal(hitPoint);
    Vector3f incomingDirection = UniformHemisphereSample(normal);
    Ray newRay = Ray(hitPoint, incomingDirection);
    Vector3f incoming = PathTrace(newRay, scene, depth + 1);
    const float p = 0.159154943092;  // 1/2pi
    float cosTheta = Dot(normal, incomingDirection);
    Vector3f BRDF = intersection.object->albedo * 0.318309886184;  // 1/pi
    return emittance + BRDF * incoming * cosTheta / p;
}

int main() {
    const int width = 500, height = 500, spp = 128;

    Vector3i **img = new Vector3i *[width];
    for (int i = 0; i < width; i++) {
        img[i] = new Vector3i[height];
    }

    Scene scene;
    scene.add(new Plane(Vector3f(0, -1, 0), Vector3f(0, 1, 0), Vector3f(1), Vector3f(0)));
    scene.add(new Sphere(Vector3f(0, 0, -1), 1, Vector3f(1), Vector3f(0)));

    clock_t startTime = clock();

    #pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < width; x++) {
        fprintf(stdout, "\rRendering... %0.2f%%", 100.0F * x / (float)width);

        for (int y = 0; y < height; y++) {
            Vector3f fragColor;

            for (int s = 0; s < spp; s++) {
                Ray viewRay;
                viewRay.origin = Vector3f(0, 0.0, 1);
                viewRay.direction = Normalize(Vector3f(2.0F * (x + random) / width - 1.0F,
                                                       1.0F - 2.0F * (y + random) / height, -1.0F));
                Vector3f color = PathTrace(viewRay, scene, 0);
                fragColor = fragColor + color / spp;
            }
            fragColor = fragColor * 1.0F;
            fragColor = Saturate(Vector3f(1) - Exp(-fragColor));
            fragColor = Pow(fragColor, 1.0F / 2.2F);
            img[x][y] = Vector3i(fragColor * 255.0F);
        }
    }

    clock_t endTime = clock();
    double time = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    printf("\rRender time: %0.3fs\n", time);

    FILE *file = fopen("Render.ppm", "w");
    fprintf(file, "P3\n%d %d\n%d\n ", width, height, 255);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            fprintf(file, "%d %d %d ", img[x][y].x, img[x][y].y, img[x][y].z);
        }
        fprintf(file, "\n");
    }
    fclose(file);

    return 0;
}