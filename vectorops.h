#pragma once

#include <math.h>
#include <vector>
#define M_PI 3.14159265359f

using namespace std;

struct Point {

    struct { double x, y, z; };
    double p[3];
    int idx = 1;

    inline Point operator*(const Point& v) const { return Point{ x * v.x, y * v.y, z * v.z }; }
    inline Point operator*(float a) const { return Point{ x * a, y * a, z * a }; }
    inline void operator+=(const Point& v) { x += v.x; y += v.y; z += v.z; }

    inline Point operator-(const Point& v) const { return Point{ x - v.x, y - v.y, z - v.z }; }

    void norm() {
        double norm = 1.0 / sqrtf(x * x + y * y + z * z);
        x *= norm; y *= norm; z *= norm;
    }

    double normd() {
        double norm = 1.0 / sqrtf(x * x + y * y + z * z);
       return norm;
    }

    Point unit() {
        double norm = 1.0f / sqrtf(x * x + y * y + z * z);
        return Point{ x * norm, y * norm, z *= norm };
    }

    Point operator+(const Point& v) const { 
        return Point{ x + v.x, y + v.y, z + v.z }; 
    }
};

Point cross(const Point& p1, const Point& p2) {
    return Point{ p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x };
}

double dot(const Point& p1, const Point& p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

double get3x3Determinant(std::vector<std::vector<double>>& m)
{
    double det = 0;
    std::vector<std::vector<double>> extendedDet(3, std::vector<double>(5));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            extendedDet[i][j] = m[i][j % 3];
        }
    }
    for (int i = 0; i < 3; i++)
    {
        double downDiagonalProduct = 1;
        double upDiagonalProduct = 1;
        for (int j = 0; j < 3; j++)
        {
            downDiagonalProduct *= extendedDet[j][j + i];
            upDiagonalProduct *= extendedDet[j][4 - j - i];
        }
        det += downDiagonalProduct - upDiagonalProduct;
    }
    return det;
}

struct f2 {
    float x, y;
};

struct Matrix4 {
    double m[16];

    Point operator *(const Point& v) {
        double x, y, z, r = 1.0f / (m[12] * v.x + m[13] * v.y + m[14] * v.z + m[15]);
        x = m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3];
        y = m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7];
        z = m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11];
        return Point{ x * r,y * r,z };
    }
};

inline
Matrix4 translate(const Point& v) {
    return Matrix4{
        1, 0, 0, v.x,
        0, 1, 0, v.y,
        0, 0, 1, v.z,
        0, 0, 0, 1 };
}

Matrix4 rotateY(double ang) {
    ang = (ang * M_PI) / 180;

    return Matrix4{
    cosf(ang),	0,	sinf(ang),	0,
    0,	1,	0,	0,
    -sinf(ang),	0,	cosf(ang),	0,
    0,	0,	0,	1 };
}

Matrix4 rotateZ(double ang) {
    ang = (ang * M_PI) / 180;

    return Matrix4{
    cosf(ang),	-sinf(ang),	0,	0,
    sinf(ang),	cosf(ang),	0,	0,
    0,	0,	1,	0,
    0,	0,	0,	1 };

}

Matrix4 rotateArbitrary(double ang, Point& u) {
    u.norm();

    ang = (ang * M_PI) / 180;
    return Matrix4{
        cosf(ang) + u.x * u.x * (1 - cosf(ang)), u.y * u.x * (1 - cosf(ang)) - u.z * sinf(ang),
               u.z * u.x * (1 - cosf(ang)) + u.y * sinf(ang),	0,

        u.y * u.x * (1 - cosf(ang)) + u.z * sinf(ang), cosf(ang) + u.y * u.y * (1 - cosf(ang)),
               u.z * u.y * (1 - cosf(ang)) - u.x * sinf(ang),	0,

        u.z * u.x * (1 - cosf(ang)) - u.y * sinf(ang), u.z * u.y * (1 - cosf(ang)) + u.x * sinf(ang),
               cosf(ang) + u.z * u.z * (1 - cosf(ang)),	0,

        0,0,0,1
    };
}

struct Camera {

    Point position;
    Point view;
    Point up;

    f2 resolution;
    f2 fov;

    void Transform(Matrix4 m) {
        view = m * view;
        up = m * up;
        position = m * position;
    }

    void Transform(Matrix4 m, Point reference) {
        Transform(translate(Point{ -reference.x, -reference.y, -reference.z }));

        Transform(m);

        Transform(translate(Point{ reference.x, reference.y, reference.z }));
    }

};