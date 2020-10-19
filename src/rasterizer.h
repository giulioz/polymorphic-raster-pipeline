#pragma once

#include <array>
#include <utility>

#include "math3d.h"

struct coord {
  int x;
  int y;

  coord() {}
  coord(int x, int y) : x(x), y(y) {}
};

Vector3 barycentric(const std::array<coord, 3> &pts, float px, float py) {
  Vector3 dest;
  Vector3 a = {(float)pts[2].x - (float)pts[0].x,
               (float)pts[1].x - (float)pts[0].x, (float)pts[0].x - px};
  Vector3 b = {(float)pts[2].y - (float)pts[0].y,
               (float)pts[1].y - (float)pts[0].y, (float)pts[0].y - py};
  Vector3 u = Vector3::cross(a, b);
  if (fabsf(u.z) < 1) {
    // triangle is degenerate, in this case return negative coordinates
    dest.x = -1;
    dest.y = 1;
    dest.z = 1;
  } else {
    dest.x = 1.f - (u.x + u.y) / u.z;
    dest.y = u.y / u.z;
    dest.z = u.x / u.z;
  }
  return dest;
}

template <unsigned long width, unsigned long height>
inline coord toScreenSpace(const Vector4 &v) {
  auto aspect = (float)height / width;
  int px = ((v.x / v.w) * aspect * width / 2) + width / 2;
  int py = height / 2 - ((v.y / v.w) * height / 2);
  return {px, py};
}

template <typename FragShader, unsigned long width, unsigned long height>
void drawTriangle(BufferWindow<width, height> &wnd,
                  const std::array<Vector4, 3> &positions,
                  const std::array<Vector4, 3> &normals,
                  std::array<double, width * height> &zBuffer,
                  FragShader &fragmentShader) {
  auto pa = toScreenSpace<width, height>(positions[0]);
  auto pb = toScreenSpace<width, height>(positions[1]);
  auto pc = toScreenSpace<width, height>(positions[2]);

  std::array<coord, 3> pts = {coord(pa.x, pa.y), coord(pb.x, pb.y),
                              coord(pc.x, pc.y)};
  coord bboxmin = {width - 1, height - 1};
  coord bboxmax = {0, 0};
  coord clamp = {width - 1, height - 1};
  for (int i = 0; i < 3; i++) {
    bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
    bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
    bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));
    bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
  }
  Vector3 P;
  for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
      Vector3 bc_screen = barycentric(pts, P.x, P.y);

      if (!(bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)) {
        P.z = 0;
        P.z += positions[0].z * bc_screen.x;
        P.z += positions[1].z * bc_screen.y;
        P.z += positions[2].z * bc_screen.z;
        Vector3 normal;
        normal.x = normals[0].x * bc_screen.x + normals[1].x * bc_screen.y +
                   normals[2].x * bc_screen.z;
        normal.y = normals[0].y * bc_screen.x + normals[1].y * bc_screen.y +
                   normals[2].y * bc_screen.z;
        normal.z = normals[0].z * bc_screen.x + normals[1].z * bc_screen.y +
                   normals[2].z * bc_screen.z;
        if (zBuffer[(int)(P.x + P.y * width)] < P.z) {
          zBuffer[(int)(P.x + P.y * width)] = P.z;
          uint32_t col;
          if (fragmentShader(&col, normal, P)) {
            wnd.setPixel(P.x, P.y, col);
          }
        }
      }
    }
  }
}
