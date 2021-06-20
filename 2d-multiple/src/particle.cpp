//
// Created by ybc on 2021/6/9.
//

#include "particle.h"


Particle::Particle(vec2 center, vec4 color, Material material) :
        position(random<Real>(center[0] - 0.08f, center[0] + 0.08f),
                 random<Real>(center[1] - 0.08f, center[1] + 0.08f)),
        velocity(0, 0),
        material(material),
        F(1,0,0,1),
        C(0),
        Jp(1.0f),
        color(color) {
}