#include <iostream>

#include "buffer.h"

#define WIDTH 800
#define HEIGHT 600

#define convCoord(x, y) (y * WIDTH + x)

inline uint32_t blend(uint32_t color1, uint32_t color2, uint8_t alpha) {
  uint32_t rb = color1 & 0xff00ff;
  uint32_t g = color1 & 0x00ff00;
  rb += ((color2 & 0xff00ff) - rb) * alpha >> 8;
  g += ((color2 & 0x00ff00) - g) * alpha >> 8;
  return (rb & 0xff00ff) | (g & 0xff00);
}

inline void vline(uint32_t *buf, int x, int y0, int y1, uint32_t color) {
  int i = convCoord(x, (int)y0);
  for (int y = y0; y < y1; y++) {
    double d = (abs(0.5 - pow((y - y0) / (float)(y1 - y0), 3)) / 2 + 0.5);
    buf[i] = blend(color, buf[i], d * 255);
    i += WIDTH;
  }
}

int main(int argc, char *argv[]) {
  auto wnd = BufferWindow<WIDTH, HEIGHT>();

  auto time = 0.0;
  while (!wnd.checkExit()) {
    wnd.startDraw();

    memset(wnd.pixels, 0, WIDTH * HEIGHT * 4);

    auto p = time / 500;
    auto scale = 100;
    auto off = HEIGHT / 2;
    for (int x = 0; x < WIDTH; x++) {
      auto y0 = (sin(p) * scale) + off;
      auto y1 = (sin(p + M_PI_2) * scale) + off;
      auto y2 = (sin(p + M_PI) * scale) + off;
      auto y3 = (sin(p + M_PI_2 * 3) * scale) + off;
      vline(wnd.pixels, x, y3, y0, 0xFFFFFF00);
      vline(wnd.pixels, x, y2, y3, 0xFFFF0000);
      vline(wnd.pixels, x, y1, y2, 0xFF00FF00);
      vline(wnd.pixels, x, y0, y1, 0xFF0000FF);

      p += 0.005 * sin(time / 500) * 3;
    }

    wnd.endDraw();

    time += wnd.getDeltaTime();
  }

  return 0;
}
