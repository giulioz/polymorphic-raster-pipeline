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

  Matrix4 viewMatrix = Matrix4::identity();

  auto time = 0.0;
  while (!wnd.checkExit()) {
    wnd.startDraw();
    wnd.clear();

    std::array<double, WIDTH *HEIGHT> zBuffer = {0};

    auto fragmentShader = [](uint32_t *col, Vector3 normal, Vector3 bary) {
      Vector3 light(0, 0, 1);
      float intensity = Vector3::dot(normal, light);
      unsigned char r = intensity * 255.0;
      unsigned char b = intensity * 255.0;
      unsigned char g = intensity * 255.0;
      unsigned char a = 0xFF;
      *col = a << 24 | r << 16 | g << 8 | b;
      return intensity > 0;
    };

    for (auto &&face : mesh.faces) {
      std::array<Vector4, 3> positions;
      for (size_t i = 0; i < 3; i++) {
        positions[i] =
            Vector4(mesh.vertexAttributes[face.indices[i]].position.x / 100,
                    mesh.vertexAttributes[face.indices[i]].position.y / 100,
                    mesh.vertexAttributes[face.indices[i]].position.z / 100, 1);
      }

      drawTriangle(wnd, mesh.vertexAttributes[face.indices[0]],
                   mesh.vertexAttributes[face.indices[1]],
                   mesh.vertexAttributes[face.indices[2]], positions, zBuffer,
                   fragmentShader);
    }

    wnd.endDraw();

    time += wnd.getDeltaTime();
  }

  return 0;
}
