#include <iostream>

#include "buffer.h"
#include "fileio.h"
#include "math3d.h"
#include "mesh.h"

#define WIDTH 800
#define HEIGHT 600

inline auto toScreenSpace(const Vector4 &v) {
  auto aspect = (float)HEIGHT / WIDTH;
  auto px = ((v.x / v.w) * aspect * WIDTH / 2) + WIDTH / 2;
  auto py = HEIGHT / 2 - ((v.y / v.w) * HEIGHT / 2);
  return std::tuple(px, py);
}

int main(int argc, char *argv[]) {
  auto wnd = BufferWindow<WIDTH, HEIGHT>();

  auto mesh = readObjFile("../dataset/teapot.obj");

  Matrix4 viewMatrix = Matrix4::identity();

  auto time = 0.0;
  while (!wnd.checkExit()) {
    wnd.startDraw();
    wnd.clear();

    for (auto &&vertex : mesh.vertexAttributes) {
      Vector4 position =
          Vector4(vertex.position.x / 100, vertex.position.y / 100,
                  vertex.position.z / 100, 1);

      auto [px, py] = toScreenSpace(position);
      if (px > 0 && py > 0 && px < WIDTH && py < HEIGHT) {
        wnd.setPixel(px, py, 0xff0000);
      }
    }

    wnd.endDraw();

    time += wnd.getDeltaTime();
  }

  return 0;
}
