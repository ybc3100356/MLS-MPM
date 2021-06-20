//
// Created by ybc on 2021/6/9.
//

#ifndef MLSMPM_ALGEBRA_H
#define MLSMPM_ALGEBRA_H

#include <glm/glm.hpp>
using Real = float;
using glm::ivec2;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::mat2;
using glm::transpose;

extern void svd(const mat2 &A, mat2 &U, mat2 &SIG, mat2 &V);

#endif //MLSMPM_ALGEBRA_H
