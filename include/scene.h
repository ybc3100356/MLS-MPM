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
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "particle.h"
#include "shader.h"
#include "camera.h"

#include <vector>

const unsigned int SCR_WIDTH = 1920;
const unsigned int SCR_HEIGHT = 1080;
using std::vector;

class Scene {
public:
    static Scene &getInstance() {
        static Scene instance;
        return instance;
    }

private:
    // shader
    unsigned int VAO;
    Shader shader;

    const size_t dim = 3;
    const size_t numParticles = 8192 * 4;
    const size_t numGrid = 128;

    const Real dt = 2.0e-4;
    const Real dx = 1.0f / numGrid;
    const Real inv_dx = 1.0f / dx;
    const Real p_vol = (dx * 0.3f) * (dx * 0.3f);
    const Real p_rho = 0.3;
    const Real p_mass = p_vol * p_rho;
    const Real E = 100;

private:
    vector<Particle> particles;
    Particle *particles_gpu;
    constexpr const static int particles_size = numParticles * sizeof(Particle);

    vector<vector<vector<vec3> > > grid_v;

    vector<vector<vector<Real> > > grid_m;

    Scene();

    void subStep();

public:
    Camera camera;

    Scene(Scene const &) = delete;

    ~Scene();

    void operator=(Scene const &) = delete;

    void processInput(GLFWwindow *window, double deltaTime);

    void processMouseMovement(float xoffset, float yoffset);

    void processMouseScroll(float yoffset);

    void update();

    void render();

    void loadShader(const GLchar *vertexPath, const GLchar *fragmentPath);
};


#endif //MLS_MPM_SCENE_H
