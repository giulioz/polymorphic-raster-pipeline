#pragma once

#include <cstring>
#include <vector>

#include "math3d.h"

using index_type = unsigned long;

struct Face {
  index_type indices[3];

  Face(index_type indices[3]) {
    std::memcpy(this->indices, indices, sizeof(index_type) * 3);
  }
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
