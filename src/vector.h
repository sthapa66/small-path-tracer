#pragma once

template <typename T>
struct Vector3 {
    T x, y, z;
    Vector3() : x(0), y(0), z(0) {}
    Vector3(T s) : x(s), y(s), z(s) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}
    template <typename U> explicit Vector3(const Vector3<U> &v) : x(v.x), y(v.y), z(v.z) {}
    Vector3<T> operator-() const {
        return Vector3<T>(-x, -y, -z);
    }
    Vector3<T> operator+(const T &s) const {
        return Vector3<T>(x + s, y + s, z + s);
    }
    Vector3<T> operator+(const Vector3<T> &v) const {
        return Vector3<T>(x + v.x, y + v.y, z + v.z);
    }
    Vector3<T> operator-(const T &s) const {
        return Vector3<T>(x - s, y - s, z - s);
    }
    Vector3<T> operator-(const Vector3<T> &v) const {
        return Vector3<T>(x - v.x, y - v.y, z - v.z);
    }
    Vector3<T> operator*(const T &s) const {
        return Vector3<T>(x * s, y * s, z * s);
    }
    Vector3<T> operator*(const Vector3<T> &v) const {
        return Vector3<T>(x * v.x, y * v.y, z * v.z);
    }
    Vector3<T> operator/(const T &s) const {
        return Vector3<T>(x / s, y / s, z / s);
    }
    Vector3<T> operator/(const Vector3<T> &v) const {
        return Vector3<T>(x / v.x, y / v.y, z / v.z);
    }
    friend std::ostream &operator<<(std::ostream &os, const Vector3<T> &v) {
        os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
        return os;
    }
};
typedef Vector3<int> Vector3i;
typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;
template <typename T> inline T
Dot(const Vector3<T> &v1, const Vector3<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
template <typename T> inline T
LengthSquared(const Vector3<T> &v) {
    return v.x*v.x + v.y*v.y + v.z*v.z;
}
template <typename T> inline T
Length(const Vector3<T> &v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
template <typename T> inline Vector3<T>
Normalize(const Vector3<T> &v) {
    return v / sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
template <typename T> inline Vector3<T>
Cross(const Vector3<T> &v1, const Vector3<T> &v2) {
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>(v1y*v2z - v1z*v2y, v1z*v2x - v1x*v2z, v1x*v2y - v1y*v2x);
}
template <typename T> inline Vector3<T>
Saturate(const Vector3<T> &v) {
    return Vector3<T>(Saturate(v.x), Saturate(v.y), Saturate(v.z));
}
template <typename T> inline Vector3<T>
Pow(const Vector3<T> &v, const T &x) {
    return Vector3<T>(pow(v.x, x), pow(v.y, x), pow(v.z, x));
}
template <typename T> inline Vector3<T>
Exp(const Vector3<T> &v) {
    return Vector3<T>(exp(v.x), exp(v.y), exp(v.z));
}
template <typename T> Vector3<T>
Max(const Vector3<T> &v, const T &x) {
    return Vector3<T>(std::max(x, v.x), std::max(x, v.y), std::max(x, v.z));
}