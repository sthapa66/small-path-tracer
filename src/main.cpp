#include <math.h>

#include <ctime>
#include <iostream>
#include <random>
#include <vector>

#include "atmosphere.h"
#include "material.h"
#include "object.h"
#include "ray.h"
#include "vector.h"


constexpr int skyWidth = 500, skyHeight = 500;


// ATMOSPHERE MODEL CONSTANTS
constexpr float gr = 6360e3;
constexpr float ar = 6420e3;
constexpr float HR = 7994;
constexpr float HM = 1200;
const Vector3f betaR(0.00000530482236632, 0.0000123952627965, 0.0000302618720618);
const Vector3f betaM(0.000021);
const Vector3f sd = Normalize(Vector3f(2, 0.5, 0.5));
const Atmosphere atmosphere(gr, ar, HR, HM, betaR, betaM, sd);


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

bool RaySphereIntersection(const Ray &ray, const float &r, float &t0, float &t1) {
    Vector3f L = -ray.origin;
    float tc = Dot(L, ray.direction);
    double d2 = Dot(L, L) - tc * tc;
    double r2 = r * r;
    if (d2 > r2) return false;
    float t1c = std::sqrt(r2 - d2);
    t0 = tc - t1c;
    t1 = tc + t1c;
    return true;
}

Vector3f UniformHemisphereSample(Vector3f n) {
    float r = sqrt(random) * 1.57079632679F;
    float t = random * 6.28318530718F;
    float s = sin(r);
    float x = s * cos(t), y = s * sin(t), z = cos(r);
    Vector3f T = Normalize(Cross(n, Vector3f(0, 1, 1)));
    Vector3f B = Cross(n, T);
    return Normalize(T * x + B * y + n * z);
}

float t0, t1;
Ray sunRay = Ray(Vector3f(0, 6360e3 + 500, 0), sd);
bool x = RaySphereIntersection(sunRay, 6420e3, t0, t1);
Vector3f sunT = ComputeTransmission(sunRay, 0.0, t1);

Vector3f PathTrace(const Ray &ray, const Scene &scene, const Vector3f const* const*skyRadianceMap, const int &depth) {
    if (depth > 4) return Vector3f(0);

    Intersection intersection = scene.intersect(ray);

    if (intersection.object == nullptr) {
        float x = ray.direction.x;
        float y = ray.direction.z;
        float z = -ray.direction.y;
        int X = (int)((skyWidth - 1) * ((x / (1.0F - z)) * 0.5F + 0.5F));
        int Y = (int)((skyHeight - 1) * ((y / (1.0F - z)) * 0.5F + 0.5F));
        if (X < 0 || X >= skyWidth || Y < 0 || Y >= skyHeight) return Vector3f(0);
        return skyRadianceMap[X][Y];
    }

    Material material = intersection.object->material;

    Vector3f hitPoint = ray.origin + ray.direction * intersection.t;
    Vector3f normal = intersection.object->getNormal(hitPoint);

    Ray shadowRay = Ray(hitPoint + normal * 0.001, sd);
    float sunVis = scene.intersect(shadowRay).object == nullptr ? 1.0f : 0.0f;
    float sunNdotL = std::max(0.0f, Dot(normal, sd));

    Vector3f incomingDirection = UniformHemisphereSample(normal);
    Ray newRay = Ray(hitPoint, incomingDirection);
    Vector3f incoming = PathTrace(newRay, scene, skyRadianceMap, depth + 1);
    const float pdf = 0.159154943092;  // 1/2pi
    float cosTheta = Dot(normal, incomingDirection);
    Vector3f BRDF = material.albedo * 0.318309886184;  // 1/pi
    return material.emittance + BRDF * sunT * sunVis * 2.0 * sunNdotL + BRDF * incoming * cosTheta / pdf;
}

void ComputeSkyRadianceMap(Vector3f **&skyRadianceMap, const int &width, const int &height) {
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < width; i++) {
        float x = i / (float)width * 2.0F - 1.0F;
        float xx = x * x;
        for (int j = 0; j < height; j++) {
            float y = j / (float)height * 2.0F - 1.0F;
            float yy = y * y;
            float denominator = 1.0F + xx + yy;
            float X = 2.0F * x / denominator;
            float Y = 2.0F * y / denominator;
            float Z = (-1.0F + xx + yy) / denominator;
            Ray viewRay = Ray(Vector3f(0, 6360e3 + 500, 0), Vector3f(X, -Z, Y));
            Vector3f skyRadiance = ComputeAtmosphereicScattering(viewRay, atmosphere);
            skyRadianceMap[i][j] = skyRadiance;
        }
    }
}

int main() {

    Vector3f **skyRadianceMap = new Vector3f*[skyWidth];
    for (int i = 0; i < skyWidth; i++) {
        skyRadianceMap[i] = new Vector3f[skyHeight];
    }
    ComputeSkyRadianceMap(skyRadianceMap, skyWidth, skyHeight);

    const int width = 500, height = 500, spp = 500;

    Vector3i **img = new Vector3i *[width];
    for (int i = 0; i < width; i++) {
        img[i] = new Vector3i[height];
    }

    Material white (Vector3f(0.8, 0.8, 0.8), 1, 0, Vector3f(0));
    Material red   (Vector3f(0.8, 0.05, 0.05), 1, 0, Vector3f(0));
    Material green (Vector3f(0.05, 0.8, 0.05), 1, 0, Vector3f(0));
    Material blue  (Vector3f(0.05, 0.05, 0.8), 1, 0, Vector3f(0));
    Material black (Vector3f(0.05, 0.05, 0.05), 1, 0, Vector3f(0));
    Material emmit (Vector3f(0.8, 0.8, 0.8), 1, 0, Vector3f(0.5));
    Material orange (Vector3f(0.05, 0.8, 0.05), 1, 0, Vector3f(0));


    Scene scene;
    scene.add(new Plane(Vector3f(0, -1, 0), Vector3f(0, 1, 0), white));
    scene.add(new Sphere(Vector3f(0, 0, -2), 1, green));

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
                Vector3f color = PathTrace(viewRay, scene, skyRadianceMap, 0);
                fragColor = fragColor + color / spp;
            }
            fragColor = fragColor * 0.7f;
            fragColor = Saturate(Vector3f(1.0) - Exp(-fragColor));
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