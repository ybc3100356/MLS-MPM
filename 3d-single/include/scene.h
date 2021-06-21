//
// Created by ybc on 2021/6/11.
//

#ifndef MLS_MPM_SCENE_H
#define MLS_MPM_SCENE_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "particle.h"
#include "shader.h"
#include "camera.h"

#include <vector>

const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 1000;
using std::vector;

class Scene {
public:
    static Scene &getInstance() {
        static Scene instance;
        return instance;
    }

    static const size_t dim = 3;
    static const size_t numParticles = 8192 * 4;
    static const size_t numGrid = 64;
    static const size_t steps = 10;

    constexpr static const Real dt = 4.0e-3;
    constexpr static const Real dx = 1.0f / numGrid;
    constexpr static const Real inv_dx = 1.0f / dx;
    constexpr static const Real p_vol = (dx * 2.0f) * (dx * 2.0f);
    constexpr static const Real p_rho = 1;
    constexpr static const Real p_mass = p_vol * p_rho;
    constexpr static const Real E = 80;

    constexpr const static int grid_v_size = numGrid * numGrid * numGrid * sizeof(vec3);
    constexpr const static int grid_m_size = numGrid * numGrid * numGrid * sizeof(Real);

private:
    Camera camera;
    // shader
    unsigned int VAO;
    Shader shader;

    vector<Particle> particles;
    Particle *particles_gpu{};
    constexpr const static int particles_size = numParticles * sizeof(Particle);

    vector<vector<vector<vec3> > > grid_v;
    vec3 *grid_v_gpu{};

    vector<vector<vector<Real> > > grid_m;
    Real *grid_m_gpu{};

    Scene();

    void subStep();

public:

    Scene(Scene const &) = delete;

    ~Scene();

    void operator=(Scene const &) = delete;

    void processInput(GLFWwindow *window, double deltaTime);

    void processMouseMovement(float xOffset, float yOffset);

    void processMouseScroll(float yOffset);

    void update();

    void render();

    void loadShader(const GLchar *vertexPath, const GLchar *fragmentPath);

private:
    // gpu related
    void gpuInit();

    void gpuFree();

    void gpuUpdate();

    // computing related
    void p2g();

    void g2p();

    void gridCompute();
};


#endif //MLS_MPM_SCENE_H
