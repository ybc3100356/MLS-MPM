//
// Created by ybc on 2021/6/9.
//

#include "particle.h"


Particle::Particle() : position(random<Real>(0.2, 0.7), random<Real>(0.3, 0.9), random<Real>(0.1, 0.9)), velocity(-10, -10, 0),
                       J(1), C(mat3()) {
}
