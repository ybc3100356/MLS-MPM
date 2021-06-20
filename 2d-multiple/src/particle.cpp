//
// Created by ybc on 2021/6/9.
//

#include "particle.h"


Particle::Particle(vec2 center, vec4 color) :
        position(random<Real>(center[0] - 0.08, center[0] + 0.08),
                 random<Real>(center[1] - 0.08, center[1] + 0.08)),
        velocity(0, 0),
        F(1),
        C(0),
        J(1),
        color(color) {
}