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

using Real = float;

class Particle {
public:
    vec2 position;
    vec2 velocity;
    mat2 C;
    Real J;

    Particle();
};

#endif //MLSMPM_PARTICLE_H
