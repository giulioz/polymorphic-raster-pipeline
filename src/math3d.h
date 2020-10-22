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

  static inline Vector4 normalize(const Vector4 &v) {
    const num factor = sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
    return Vector4{v.x / factor, v.y / factor, v.z / factor, v.w / factor};
  }

  const num length() const { return sqrt(x * x + y * y + z * z + w * w); }
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

  static auto perspective(const num near, const num far) {
    return Matrix4{
        2 * near, 0.0, 0.0, 0.0,  //
        0.0, 2 * near, 0.0, 0.0,  //
        0.0, 0.0, -(far + near) / (far - near), -(2 * far * near) / (far - near),  //
        0.0, 0.0, -1.0, 0.0,  //
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

  static auto transformVector4(const Vector4 &v, const Matrix4 &m) {
    return Vector4{v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + m.m41 * v.w,
                   v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + m.m42 * v.w,
                   v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + m.m43 * v.w,
                   v.x * m.m14 + v.y * m.m24 + v.z * m.m34 + m.m44 * v.w};
  }

  Matrix4 operator*(const Matrix4 &value2) const {
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

  static Matrix4 invert(const Matrix4 &matrix) {
    Matrix4 result;

    float a = matrix.m11, b = matrix.m12, c = matrix.m13, d = matrix.m14;
    float e = matrix.m21, f = matrix.m22, g = matrix.m23, h = matrix.m24;
    float i = matrix.m31, j = matrix.m32, k = matrix.m33, l = matrix.m34;
    float m = matrix.m41, n = matrix.m42, o = matrix.m43, p = matrix.m44;

    float kp_lo = k * p - l * o;
    float jp_ln = j * p - l * n;
    float jo_kn = j * o - k * n;
    float ip_lm = i * p - l * m;
    float io_km = i * o - k * m;
    float in_jm = i * n - j * m;

    float a11 = +(f * kp_lo - g * jp_ln + h * jo_kn);
    float a12 = -(e * kp_lo - g * ip_lm + h * io_km);
    float a13 = +(e * jp_ln - f * ip_lm + h * in_jm);
    float a14 = -(e * jo_kn - f * io_km + g * in_jm);

    float det = a * a11 + b * a12 + c * a13 + d * a14;

    if (fabs(det) < FLT_EPSILON) {
      return Matrix4(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                     NAN, NAN, NAN, NAN);
    }

    float invDet = 1.0f / det;

    result.m11 = a11 * invDet;
    result.m21 = a12 * invDet;
    result.m31 = a13 * invDet;
    result.m41 = a14 * invDet;

    result.m12 = -(b * kp_lo - c * jp_ln + d * jo_kn) * invDet;
    result.m22 = +(a * kp_lo - c * ip_lm + d * io_km) * invDet;
    result.m32 = -(a * jp_ln - b * ip_lm + d * in_jm) * invDet;
    result.m42 = +(a * jo_kn - b * io_km + c * in_jm) * invDet;

    float gp_ho = g * p - h * o;
    float fp_hn = f * p - h * n;
    float fo_gn = f * o - g * n;
    float ep_hm = e * p - h * m;
    float eo_gm = e * o - g * m;
    float en_fm = e * n - f * m;

    result.m13 = +(b * gp_ho - c * fp_hn + d * fo_gn) * invDet;
    result.m23 = -(a * gp_ho - c * ep_hm + d * eo_gm) * invDet;
    result.m33 = +(a * fp_hn - b * ep_hm + d * en_fm) * invDet;
    result.m43 = -(a * fo_gn - b * eo_gm + c * en_fm) * invDet;

    float gl_hk = g * l - h * k;
    float fl_hj = f * l - h * j;
    float fk_gj = f * k - g * j;
    float el_hi = e * l - h * i;
    float ek_gi = e * k - g * i;
    float ej_fi = e * j - f * i;

    result.m14 = -(b * gl_hk - c * fl_hj + d * fk_gj) * invDet;
    result.m24 = +(a * gl_hk - c * el_hi + d * ek_gi) * invDet;
    result.m34 = -(a * fl_hj - b * el_hi + d * ej_fi) * invDet;
    result.m44 = +(a * fk_gj - b * ek_gi + c * ej_fi) * invDet;

    return result;
  }

  static Matrix4 transpose(const Matrix4 &matrix) {
    Matrix4 result;

    result.m11 = matrix.m11;
    result.m12 = matrix.m21;
    result.m13 = matrix.m31;
    result.m14 = matrix.m41;
    result.m21 = matrix.m12;
    result.m22 = matrix.m22;
    result.m23 = matrix.m32;
    result.m24 = matrix.m42;
    result.m31 = matrix.m13;
    result.m32 = matrix.m23;
    result.m33 = matrix.m33;
    result.m34 = matrix.m43;
    result.m41 = matrix.m14;
    result.m42 = matrix.m24;
    result.m43 = matrix.m34;
    result.m44 = matrix.m44;

    return result;
  }
};
