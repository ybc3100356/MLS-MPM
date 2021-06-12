//
// Created by ybc on 2021/6/11.
//

#ifndef MLS_MPM_SCENE_H
#define MLS_MPM_SCENE_H


#include "particle.h"
#include "shader.h"

#include <vector>

using std::vector;

class Scene {
    const size_t dim = 2;
    const size_t numParticles = 4096;
    const size_t numGrid = 128;

    const Real dx = 1.0f / numGrid;
    const Real inv_dx = 1.0f / dx;
    const Real p_vol = (dx * 0.5f) * (dx * 0.5f);
    const Real p_rho = 0.35;
    const Real p_mass = p_vol * p_rho;
    const Real E = 200;

    vector<Particle> particles;

    vector<vector<vec2>> grid_v;

    vector<vector<Real>> grid_m;

    void subStep();

    // shader
    unsigned int VAO;
    Shader &shader;
public:
    const Real dt = 2.0e-3;

    explicit Scene(const Shader &shader);

    ~Scene();

    void update();

    void render();
};


#endif //MLS_MPM_SCENE_H
