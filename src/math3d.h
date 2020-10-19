#pragma once

#include <algorithm>
#include <array>
#include <cmath>

using num = float;

num radians(num degrees) { return degrees / 180 * M_PI; }

struct Vector2 {
  num x, y;

  Vector2() {}

  Vector2(num x, num y) : x(x), y(y) {}

  Vector2(const Vector2 &other) : x(other.x), y(other.y) {}

  static inline Vector2 normalize(const Vector2 &v) {
    const num factor = sqrt(v.x * v.x + v.y * v.y);
    return Vector2{v.x / factor, v.y / factor};
  }

  const num length() { return sqrt(x * x + y * y); }

  const Vector2 operator+(const Vector2 &value2) const {
    return Vector2{x + value2.x, y + value2.y};
  }

  const Vector2 operator-(const Vector2 &value2) const {
    return Vector2{x - value2.x, y - value2.y};
  }

  const Vector2 operator*(num value2) const {
    return Vector2{x * value2, y * value2};
  }
};

struct Vector3 {
  num x, y, z;

  Vector3() : x(0), y(0), z(0) {}
  Vector3(num x, num y, num z) : x(x), y(y), z(z) {}
  Vector3(const Vector3 &other) : x(other.x), y(other.y), z(other.z) {}

  static inline Vector3 cross(const Vector3 &a, const Vector3 &b) {
    Vector3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
  }

  static inline num dot(const Vector3 &a, const Vector3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  static inline Vector3 reflect(const Vector3 &N, const Vector3 &I) {
    const num d = dot(N, I);
    return I - N * d * 2.0f;
  }

  static inline Vector3 max(const Vector3 &v, num a) {
    return Vector3{fmax(v.x, a), fmax(v.y, a), fmax(v.z, a)};
  }

  static inline Vector3 abs(const Vector3 &v) {
    return Vector3{fabs(v.x), fabs(v.y), fabs(v.z)};
  }

  static inline Vector3 normalize(const Vector3 &v) {
    const num factor = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return Vector3{v.x / factor, v.y / factor, v.z / factor};
  }

  const num length() const { return sqrt(x * x + y * y + z * z); }

  const Vector3 operator+(const Vector3 &value2) const {
    return Vector3{x + value2.x, y + value2.y, z + value2.z};
  }

  const Vector3 operator-(const Vector3 &value2) const {
    return Vector3{x - value2.x, y - value2.y, z - value2.z};
  }

  const Vector3 operator-() const { return Vector3{-x, -y, -z}; }

  const Vector3 operator*(const Vector3 &value2) const {
    return Vector3(x * value2.x, y * value2.y, z * value2.z);
  }

  const Vector3 operator*(num value2) const {
    return Vector3{x * value2, y * value2, z * value2};
  }

  const Vector3 operator/(num value2) const {
    return Vector3{x / value2, y / value2, z / value2};
  }
};

struct Vector4 {
  num x, y, z, w;

  Vector4() {}

  Vector4(const Vector4 &other)
      : x(other.x), y(other.y), z(other.z), w(other.w) {}

  Vector4(num _x, num _y, num _z, num _w) {
    x = _x;
    y = _y;
    z = _z;
    w = _w;
  }

  Vector4(const Vector3 &v3, num _w) {
    x = v3.x;
    y = v3.y;
    z = v3.z;
    w = _w;
  }
};

struct Matrix4 {
  num m11, m21, m31, m41;
  num m12, m22, m32, m42;
  num m13, m23, m33, m43;
  num m14, m24, m34, m44;

  Matrix4()
      : m11(0),
        m21(0),
        m31(0),
        m41(0),
        m12(0),
        m22(0),
        m32(0),
        m42(0),
        m13(0),
        m23(0),
        m33(0),
        m43(0),
        m14(0),
        m24(0),
        m34(0),
        m44(0) {}

  Matrix4(num m11, num m21, num m31, num m41, num m12, num m22, num m32,
          num m42, num m13, num m23, num m33, num m43, num m14, num m24,
          num m34, num m44)
      : m11(m11),
        m21(m21),
        m31(m31),
        m41(m41),
        m12(m12),
        m22(m22),
        m32(m32),
        m42(m42),
        m13(m13),
        m23(m23),
        m33(m33),
        m43(m43),
        m14(m14),
        m24(m24),
        m34(m34),
        m44(m44) {}

  static auto identity() {
    return Matrix4{
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    };
  }

  static auto scale(const Vector3 &scale) {
    return Matrix4{
        scale.x, 0.0, 0.0,     0.0, 0.0, scale.y, 0.0, 0.0,
        0.0,     0.0, scale.z, 0.0, 0.0, 0.0,     0.0, 1.0,
    };
  }

  static auto translation(const Vector3 &t) {
    return Matrix4{
        1.0, 0.0, 0.0, t.x, 0.0, 1.0, 0.0, t.y,
        0.0, 0.0, 1.0, t.z, 0.0, 0.0, 0.0, 1.0,
    };
  }

  static auto rotationX(const num radians) {
    auto c = cos(radians);
    auto s = sin(radians);
    return Matrix4{
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, c,    s,    0.0f,
        0.0f, -s,   c,    0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    };
  }

  static auto rotationY(const num radians) {
    auto c = cos(radians);
    auto s = sin(radians);
    return Matrix4{
        c, 0.0f, -s, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        s, 0.0f, c,  0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    };
  }

  static auto rotationZ(const num radians) {
    auto c = cos(radians);
    auto s = sin(radians);
    return Matrix4{
        c,    s,    0.0f, 0.0f, -s,   c,    0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    };
  }

  static auto perspective() {
    return Matrix4{
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    };
  }

  static auto transformVector3(const Vector3 &v, const Matrix4 &m,
                               num w = 1.0) {
    return Vector3{v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + m.m41 * w,
                   v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + m.m42 * w,
                   v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + m.m43 * w};
  }

  static auto transformVector4(const Vector3 &v, const Matrix4 &m,
                               num w = 1.0) {
    return Vector4{v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + m.m41 * w,
                   v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + m.m42 * w,
                   v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + m.m43 * w,
                   v.x * m.m14 + v.y * m.m24 + v.z * m.m34 + m.m44 * w};
  }

  Matrix4 operator*(Matrix4 value2) {
    Matrix4 m;

    m.m11 = m11 * value2.m11 + m12 * value2.m21 + m13 * value2.m31 +
            m14 * value2.m41;
    m.m12 = m11 * value2.m12 + m12 * value2.m22 + m13 * value2.m32 +
            m14 * value2.m42;
    m.m13 = m11 * value2.m13 + m12 * value2.m23 + m13 * value2.m33 +
            m14 * value2.m43;
    m.m14 = m11 * value2.m14 + m12 * value2.m24 + m13 * value2.m34 +
            m14 * value2.m44;

    m.m21 = m21 * value2.m11 + m22 * value2.m21 + m23 * value2.m31 +
            m24 * value2.m41;
    m.m22 = m21 * value2.m12 + m22 * value2.m22 + m23 * value2.m32 +
            m24 * value2.m42;
    m.m23 = m21 * value2.m13 + m22 * value2.m23 + m23 * value2.m33 +
            m24 * value2.m43;
    m.m24 = m21 * value2.m14 + m22 * value2.m24 + m23 * value2.m34 +
            m24 * value2.m44;

    m.m31 = m31 * value2.m11 + m32 * value2.m21 + m33 * value2.m31 +
            m34 * value2.m41;
    m.m32 = m31 * value2.m12 + m32 * value2.m22 + m33 * value2.m32 +
            m34 * value2.m42;
    m.m33 = m31 * value2.m13 + m32 * value2.m23 + m33 * value2.m33 +
            m34 * value2.m43;
    m.m34 = m31 * value2.m14 + m32 * value2.m24 + m33 * value2.m34 +
            m34 * value2.m44;

    m.m41 = m41 * value2.m11 + m42 * value2.m21 + m43 * value2.m31 +
            m44 * value2.m41;
    m.m42 = m41 * value2.m12 + m42 * value2.m22 + m43 * value2.m32 +
            m44 * value2.m42;
    m.m43 = m41 * value2.m13 + m42 * value2.m23 + m43 * value2.m33 +
            m44 * value2.m43;
    m.m44 = m41 * value2.m14 + m42 * value2.m24 + m43 * value2.m34 +
            m44 * value2.m44;

    return m;
  }
};
