#pragma once

#include <SDL.h>

inline uint32_t alphaBlend(uint32_t color1, uint32_t color2, uint8_t alpha) {
  uint32_t rb = color1 & 0xff00ff;
  uint32_t g = color1 & 0x00ff00;
  rb += ((color2 & 0xff00ff) - rb) * alpha >> 8;
  g += ((color2 & 0x00ff00) - g) * alpha >> 8;
  return (rb & 0xff00ff) | (g & 0xff00);
}

template <size_t width, size_t height>
class BufferWindow {
 private:
  SDL_Window *window;
  SDL_Surface *screenSurface;

 public:
  using color_type = uint32_t;

  color_type *pixels;
  uint64_t last = SDL_GetPerformanceCounter();

  BufferWindow() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
      printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
    } else {
      window = SDL_CreateWindow("SDL Output", SDL_WINDOWPOS_UNDEFINED,
                                SDL_WINDOWPOS_UNDEFINED, width, height,
                                SDL_WINDOW_SHOWN);
      if (window == NULL) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
      } else {
        screenSurface = SDL_GetWindowSurface(window);
        pixels = (color_type *)screenSurface->pixels;
      }
    }
  }

  ~BufferWindow() {
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  inline size_t convCoord(size_t x, size_t y) { return y * width + x; }

  inline void setPixel(size_t i, color_type color) { pixels[i] = color; }

  inline void setPixel(size_t x, size_t y, color_type color) {
    pixels[convCoord(x, y)] = color;
  }

  void clear() { memset(pixels, 0, width * height * sizeof(color_type)); }

  auto checkExit() {
    SDL_Event e;
    while (SDL_PollEvent(&e) != 0) {
      if (e.type == SDL_QUIT) {
        return true;
      }
    }

    return false;
  }

  auto getDeltaTime() {
    auto now = SDL_GetPerformanceCounter();
    auto dt =
        (double)((now - last) * 1000 / (double)SDL_GetPerformanceFrequency());
    last = now;
    return dt;
  }

  void startDraw() { SDL_LockSurface(screenSurface); }

  void endDraw() {
    SDL_UnlockSurface(screenSurface);
    SDL_UpdateWindowSurface(window);
  }
};
