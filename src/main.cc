#include <array>
#include <iostream>
#include <utility>

#include "buffer.h"
#include "fileio.h"
#include "math3d.h"
#include "mesh.h"
#include "rasterizer.h"

#define WIDTH 800
#define HEIGHT 600

int main(int argc, char *argv[]) {
  auto wnd = BufferWindow<WIDTH, HEIGHT>();

  auto mesh = readObjFile("../dataset/teapot.obj");

  std::array<double, WIDTH *HEIGHT> zBuffer = {0};

  auto time = 0.0;
  while (!wnd.checkExit()) {
    wnd.startDraw();
    wnd.clear();
    std::fill(zBuffer.begin(), zBuffer.end(), 0);

    auto fragmentShader = [](uint32_t *col, Vector3 normal, Vector3 bary) {
      Vector3 light(0, 0, -1);
      float intensity = Vector3::dot(normal, light);
      unsigned char r = intensity * 255.0;
      unsigned char b = intensity * 255.0;
      unsigned char g = intensity * 255.0;
      unsigned char a = 0xFF;
      *col = a << 24 | r << 16 | g << 8 | b;
      return intensity > 0;
    };

    Matrix4 transform =
        Matrix4::rotationY(time / 1000.0) * Matrix4::rotationX(time / 1000.0) *
        Matrix4::scale(Vector3(1.0 / 100, 1.0 / 100, 1.0 / 100)) *
        Matrix4::translation(Vector3(0, 0, 1.5));

    Matrix4 transformInverse = Matrix4::transpose(Matrix4::invert(transform));

    Matrix4 projection = Matrix4::perspective(1, 5);

    for (auto &&face : mesh.faces) {
      std::array<Vector4, 3> positions;
      std::array<Vector4, 3> normals;
      for (size_t i = 0; i < 3; i++) {
        positions[i] = Matrix4::transformVector4(
            mesh.vertexAttributes[face.indices[i]].position, transform);
        positions[i] = Matrix4::transformVector4(positions[i], projection);
        positions[i].x /= positions[i].w;
        positions[i].y /= positions[i].w;
        positions[i].z /= positions[i].w;

        normals[i] = Vector4::normalize(Matrix4::transformVector4(
            mesh.vertexAttributes[face.indices[i]].normal, transformInverse));
      }

      drawTriangle(wnd, positions, normals, zBuffer, fragmentShader);
    }

    wnd.endDraw();

    time += wnd.getDeltaTime();
  }

  return 0;
}
