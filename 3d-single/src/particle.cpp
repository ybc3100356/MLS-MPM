//
// Created by ybc on 2021/6/9.
//

#include "particle.h"


Particle::Particle() : position(
        random<Real>(0.05, 0.55),
        random<Real>(0.05, 0.55),
        random<Real>(0.05, 0.55)),
        velocity(0, -1, 0),
                       J(1), C(mat3(0)) {
}
