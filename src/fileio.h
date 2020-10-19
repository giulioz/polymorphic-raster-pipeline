#pragma once

#include <array>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "mesh.h"

Mesh<TexCoordNormalVertex> readObjFile(const std::string& filename) {
  Mesh<TexCoordNormalVertex> mesh;
  std::vector<Vector3> pos;
  std::vector<Vector3> norm;
  std::vector<Vector3> tcoord;

  std::ifstream infile(filename);
  std::string line;
  while (std::getline(infile, line)) {
    if (line.rfind("v ", 0) == 0) {
      std::string substr = line.substr(1);
      std::stringstream ss(substr);
      double vals[3];
      ss >> vals[0];
      ss >> vals[1];
      ss >> vals[2];
      pos.push_back(Vector3(vals[0], vals[1], vals[2]));
    } else if (line.rfind("vn ", 0) == 0) {
      std::string substr = line.substr(2);
      std::stringstream ss(substr);
      double vals[3];
      ss >> vals[0];
      ss >> vals[1];
      ss >> vals[2];
      norm.push_back(Vector3(vals[0], vals[1], vals[2]));
    } else if (line.rfind("vt ", 0) == 0) {
      std::string substr = line.substr(2);
      std::stringstream ss(substr);
      double vals[3];
      ss >> vals[0];
      ss >> vals[1];
      ss >> vals[2];
      tcoord.push_back(Vector3(vals[0], vals[1], vals[2]));
    } else if (line.rfind("f ", 0) == 0) {
      std::array<unsigned long, 3> vals;
      std::string substr = line.substr(1);
      std::stringstream ss(substr);
      int i = 0;
      while (ss.peek() == ' ' || ss.peek() == '/') ss.ignore();
      while (ss >> vals[i]) {
        vals[i]--;
        i++;
        while (ss.peek() == ' ' || ss.peek() == '/') ss.ignore();
      }
      mesh.faces.push_back(Face(vals));
    }
  }

  for (size_t i = 0; i < pos.size(); i++) {
    TexCoordNormalVertex v;
    v.position = pos[i];
    if (norm.size() == pos.size()) {
      v.normal = norm[i];
    }
    if (tcoord.size() == pos.size()) {
      v.tcoord = tcoord[i];
    }
    mesh.vertexAttributes.push_back(v);
  }

  return mesh;
}
