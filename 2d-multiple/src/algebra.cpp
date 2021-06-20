//
// Created by ybc on 2021/6/9.
//

#include "algebra.h"

inline Real determinant(const mat2 &mat) {
    return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

#include <cmath>
#include <vector>

using std::sqrt;
using std::atan2;
using std::cos;
using std::sin;
using std::fabs;
using std::hypot;
using std::vector;

void svd(const mat2 &A, mat2 &U, mat2 &SIG, mat2 &V) {
    mat2 R, S;
    Real x = A[0][0] + A[1][1], y = A[0][1] - A[1][0];
    Real scale = (1.0f / sqrt(x * x + y * y));
    Real cosx = x * scale;
    Real sinx = y * scale;
    R = mat2(cosx, sinx, -sinx, cosx);
    S = transpose(R) * A;
    Real c = 0, s = 0;
    Real s1 = 0, s2 = 0;
    if (fabs(S[1][0]) < 1e-5f) {
        c = 1;
        s = 0;
        s1 = S[0][0];
        s2 = S[1][1];
    } else {
        Real tao = 0.5f * (S[0][0] - S[1][1]);
        Real w = sqrt(tao * tao + S[1][0] * S[1][0]);
        Real t = 0;
        if (tao > 0) {
            t = S[1][0] / (tao + w);
        } else {
            t = S[1][0] / (tao - w);
        }
        c = 1.0f / sqrt(t * t + 1);
        s = -t * c;
        s1 = c * c * S[0][0] - 2 * c * s * S[1][0] + s * s * S[1][1];
        s2 = s * s * S[0][0] + 2 * c * s * S[1][0] + c * c * S[1][1];
    }
    if (s1 < s2) {
        Real tmp = s1;
        s1 = s2;
        s2 = tmp;
        V = mat2(-s, -c, c, -s);
    } else {
        V = mat2(c, -s, s, c);
    }
    U = R * V;
    SIG = mat2(s1, 0, 0, s2);
}