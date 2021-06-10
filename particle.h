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

    Particle(vec2<Real> &position, vec2<Real> &velocity) : position(position), velocity(velocity) {}

    Particle() {
        position[0] = random<Real>(0.2, 0.8);
        position[1] = random<Real>(0.2, 0.8);
        velocity[0] = random<Real>(-0.8, 0.8);
        velocity[1] = random<Real>(-0.8, 0.8);
//        velocity[0] = 1.0;
//        velocity[1] = 0.0;
    }
};


#endif //MLSMPM_PARTICLE_H
