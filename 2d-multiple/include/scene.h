//
// Created by ybc on 2021/6/11.
//

#ifndef MLS_MPM_SCENE_H
#define MLS_MPM_SCENE_H

#ifdef __CUDACC__
#define CUDA_HOST_DEV __global__
#else
#define CUDA_HOST_DEV
#endif

#include "particle.h"
#include "shader.h"

#include <vector>

using std::vector;

class Scene {
public:
    const static size_t dim = 2;
    const static size_t numParticlesPerObject = 8192 ;
    const static size_t numObject = 3;
    const static size_t numParticles = numParticlesPerObject * numObject;
    const static size_t numGrid = 256;

    constexpr const static Real E = 5e3;
    constexpr const static Real nu = 0.2;
    constexpr const static Real mu_0 = E / (2.0f * (1.0f + nu));
    constexpr const static Real lambda_0 = E * nu / ((1.0f + nu) * (1.0f - 2.0f * nu)); // Lame parameters
    constexpr const static Real dt = 2.0e-4;
    constexpr const static size_t steps = 25;

    constexpr const static Real dx = 1.0f / numGrid;
    constexpr const static Real inv_dx = 1.0f / dx;
    constexpr const static Real p_vol = (dx * 0.5f) * (dx * 0.5f);
    constexpr const static Real p_rho = 1;
    constexpr const static Real p_mass = p_vol * p_rho;

    constexpr const static int grid_v_size = numGrid * numGrid * sizeof(vec2);
    constexpr const static int grid_m_size = numGrid * numGrid * sizeof(Real);

private:
    vector<Particle> particles;
    Particle *particles_gpu{};
    constexpr const static int particles_size = numParticles * sizeof(Particle);

    vector<vector<vec2>> grid_v;
    vec2 *grid_v_gpu{};

    vector<vector<Real>> grid_m;
    Real *grid_m_gpu{};

    void subStep();

    // shader
    unsigned int VAO;
    Shader &shader;
public:

    explicit Scene(const Shader &shader);

    ~Scene();

    void update();

    void render();

    void gpuInit();

    void gpuFree();

    void gpuUpdate();

    void p2g();

    void g2p();

    void gridCompute();
};


#endif //MLS_MPM_SCENE_H
