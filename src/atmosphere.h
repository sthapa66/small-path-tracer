#pragma once

#include "ray.h"
#include "vector.h"

struct Atmosphere {
    float gr, ar, HR, HM;
    Vector3f betaR, betaM, sunDirection;
    Atmosphere(float GR, float AR, float hr, float hm, Vector3f BetaR, Vector3f BetaM,
               Vector3f SunDirection) {
        gr = GR;
        ar = AR;
        HR = hr;
        HM = hm;
        betaR = BetaR;
        betaM = BetaM;
        sunDirection = SunDirection;
    }
};

bool RaySphereIntersection(const Vector3f &ro, const Vector3f &rd, const float &r, float &t0,
                           float &t1) {
    Vector3f L = -ro;
    float tc = Dot(L, rd);
    double d2 = Dot(L, L) - tc * tc;
    double r2 = r * r;
    if (d2 > r2) return false;
    float t1c = std::sqrt(r2 - d2);
    t0 = tc - t1c;
    t1 = tc + t1c;
    return true;
}

Vector3f ComputeTransmission(const Ray &ray, const float &t0, const float &t1) {
    Vector3f betaR(0.00000530482236632, 0.0000123952627965, 0.0000302618720618);
    Vector3f betaM(0.000021);
    float opticalDepthR = 0;
    float opticalDepthM = 0;
    float ds = (t1 - t0) / 8.0F;
    for (int i = 0; i < 8; i++) {
        Vector3f samplePosition = ray.origin + ray.direction * ds * ((float)i + 0.5F);
        float height = Length(samplePosition) - 6360e3;
        opticalDepthR += exp(-height / 7994.0F) * ds;
        opticalDepthM += exp(-height / 1200.0F) * ds;
    }
    return Exp(-betaR * opticalDepthR - betaM * opticalDepthM);
}

Vector3f ComputeAtmosphereicScattering(Ray ray, Atmosphere atmosphere) {
    float t0 = 0, t1 = 0;
    if (!RaySphereIntersection(ray.origin, ray.direction, 6420e3, t0, t1) || t1 < 0.0F) {
        return Vector3f(1, 0, 0);
    }
    t0 = t0 > 0.0F ? t0 : 0.0F;

    // rayleigh phase function
    float mu = Dot(ray.direction, atmosphere.sunDirection);
    float mu2 = mu * mu;
    float phaseR = 0.0596831036595F * (1.0F + mu2);

    // mie phase function
    float g = 0.76;
    float g2 = g * g;
    float phaseM =
        0.119366207319 * ((1 - g2) * (1 + mu2)) / ((2 + g2) * pow(1 + g2 - 2 * g * mu, 1.5));

    float opticalDepthR = 0;
    float opticalDepthM = 0;

    Vector3f skyColor = Vector3f(0);

    int nSamples = 32;
    float segmentLength = (t1 - t0) / float(nSamples);

    for (int i = 0; i < nSamples; i++) {
        Vector3f samplePosition =
            ray.origin + ray.direction * t0 + ray.direction * segmentLength * (float(i) + 0.5F);
        float height = Length(samplePosition) - atmosphere.gr;
        float hr = exp(-height / atmosphere.HR) * segmentLength;
        float hm = exp(-height / atmosphere.HM) * segmentLength;
        opticalDepthR += hr;
        opticalDepthM += hm;

        float t0l, t1l;
        RaySphereIntersection(samplePosition, atmosphere.sunDirection, atmosphere.ar, t0l, t1l);

        int nSamplesl = 8;
        float segmentLengthl = t1l / float(nSamplesl);
        float opticalDepthRl = 0;
        float opticalDepthMl = 0;
        for (int j = 0; j < nSamples; j++) {
            Vector3f samplePositionl =
                samplePosition + atmosphere.sunDirection * segmentLengthl * (float(j) + 0.5F);
            float heightl = Length(samplePositionl) - atmosphere.gr;
            opticalDepthRl += exp(-heightl / atmosphere.HR) * segmentLengthl;
            opticalDepthMl += exp(-heightl / atmosphere.HM) * segmentLengthl;
        }
        Vector3f color =
            Exp(-(atmosphere.betaR * opticalDepthR + atmosphere.betaM * opticalDepthM)) *
            Exp(-(atmosphere.betaR * opticalDepthRl + atmosphere.betaM * opticalDepthMl)) *
            (atmosphere.betaR * hr * phaseR + atmosphere.betaM * hm * phaseM);
        skyColor = skyColor + color;
    }
    return skyColor * 30;
}
