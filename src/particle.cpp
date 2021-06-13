//
// Created by ybc on 2021/6/9.
//

#include "particle.h"


Particle::Particle() : position(random<Real>(0.2, 0.7), random<Real>(0.3, 0.9)), velocity(0,-10),
                       J(1), C() {
}
