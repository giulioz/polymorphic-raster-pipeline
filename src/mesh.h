#pragma once

#include <array>
#include <cstring>
#include <vector>

#include "math3d.h"

using index_type = unsigned long;

struct Face {
  std::array<index_type, 3> indices;

  template <typename T>
  Face(const T &indices) : indices(indices) {}
};

struct TexCoordNormalVertex {
  Vector3 position;
  Vector3 normal;
  Vector3 tcoord;
};

template <typename VertexT = TexCoordNormalVertex>
struct Mesh {
  std::vector<VertexT> vertexAttributes;
  std::vector<Face> faces;

  Mesh() {}

  Mesh(const std::vector<VertexT> &vertexAttributes,
       const std::vector<Face> &faces)
      : vertexAttributes(vertexAttributes), faces(faces) {}

  Mesh(const Mesh<VertexT> &copy)
      : vertexAttributes(copy.vertexAttributes), faces(copy.faces) {}
};
