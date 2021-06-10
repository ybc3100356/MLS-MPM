//
// Created by ybc on 2021/6/9.
//

#ifndef MLSMPM_PARTICLE_H
#define MLSMPM_PARTICLE_H

#include "algebra.h"
#include <random>

template<typename Numeric, typename Generator = std::mt19937>
Numeric random(Numeric from, Numeric to) {
    thread_local static Generator gen(std::random_device{}());

    using dist_type = typename std::conditional
            <std::is_integral<Numeric>::value, std::uniform_int_distribution<Numeric>, std::uniform_real_distribution<Numeric> >::type;

    thread_local static dist_type dist;

    return dist(gen, typename dist_type::param_type{from, to});
}

using Real = double;

class Particle {
public:
    vec2<Real> position;
    vec2<Real> velocity;
    mat2<Real> C;
    Real J;

    Particle(vec2<Real> &position, vec2<Real> &velocity, mat2<Real> &C)
            : position(position), velocity(velocity), J(1.0) {
        C = mat2<Real>();
    }

    Particle() {
        position[0] = random<Real>(0.1, 0.6); // x
        position[1] = random<Real>(0.4, 0.7); // y
        velocity = vec2<Real>{-3, -3};
        J = 1;
        C = mat2<Real>();
    }
};


#endif //MLSMPM_PARTICLE_H
