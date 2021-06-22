//
// Created by ybc on 2021/6/11.
//

#include "scene.h"
#include <cassert>
#include <cooperative_groups.h>

const int numBlock = 128;
const int numThread = 128;

__global__ void gpuCompute(Particle *particles, vec3 *grid_v, Real *grid_m) {
    cooperative_groups::grid_group grid = cooperative_groups::this_grid();
    int threadId = int(grid.thread_rank());
    int totalThreadNum = int(grid.size());

    for (int step = 0; step < Scene::steps; step++) {
        // has to memset on each iter start
        if (threadId == 0) {
            memset(grid_v, 0, Scene::grid_v_size);
            memset(grid_m, 0, Scene::grid_m_size);
        }
        grid.sync();

        // p2g
//        assert(Scene::numParticles % totalThreadNum == 0);
        int pgRepeatTime = int(Scene::numParticles / totalThreadNum);
        for (int i = 0; i < pgRepeatTime; i++) {
            size_t idx = threadId * pgRepeatTime + i;
            if (idx > Scene::numParticles) continue;
            auto base = (ivec3) (particles[idx].position * Scene::inv_dx - 0.5f);
            auto fx = (particles[idx].position * Scene::inv_dx) - (vec3) base;
            // quadratic B-spline weights
            vec3 w[] = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                        0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                        0.5f * (fx - 0.5f) * (fx - 0.5f)};
            Real stress =
                    -Scene::dt * Scene::p_vol * (particles[idx].J - 1) * 4 * Scene::inv_dx * Scene::inv_dx * Scene::E;
            mat3 affine =
                    mat3(vec3(stress, 0, 0), vec3(0, stress, 0), vec3(0, 0, stress)) + Scene::p_mass * particles[idx].C;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        auto offset = ivec3(i, j, k);
                        auto dpos = ((vec3) offset - fx) * Scene::dx;
                        auto weight = w[i][0] * w[j][1] * w[k][2];
                        auto index = base + offset;
                        if (!(index[0] < Scene::numGrid && index[1] < Scene::numGrid &&
                              index[2] < Scene::numGrid))
                            continue;
                        auto dv = weight * (Scene::p_mass * particles[idx].velocity + affine * dpos);
                        auto target_idx = index[0] * Scene::numGrid * Scene::numGrid + index[1] * Scene::numGrid + index[2];
                        atomicAdd(&(grid_v[target_idx][0]), dv[0]);
                        atomicAdd(&(grid_v[target_idx][1]), dv[1]);
                        atomicAdd(&(grid_v[target_idx][2]), dv[2]);
                        atomicAdd(&(grid_m[target_idx]), weight * Scene::p_mass);
                    }
                }
            }
        }
        grid.sync();

        // grid
        assert((Scene::numGrid * Scene::numGrid * Scene::numGrid) % totalThreadNum == 0);
        int gridRepeatTime = int((Scene::numGrid * Scene::numGrid * Scene::numGrid) / totalThreadNum);
        for (int i = 0; i < gridRepeatTime; i++) {
            size_t idx = threadId * gridRepeatTime + i;
            if (grid_m[idx] > 0) {
                auto inv_m = 1.0f / grid_m[idx];
                grid_v[idx] = inv_m * grid_v[idx];
                grid_v[idx][1] -= Scene::dt * 9.8;
                auto bound = 3;

                size_t i = idx / (Scene::numGrid * Scene::numGrid);
                size_t j = (idx / Scene::numGrid) % Scene::numGrid;
                size_t k = idx % Scene::numGrid;

                if (i < bound && grid_v[idx].x < 0)
                    grid_v[idx].x = 0;
                if (i > Scene::numGrid - bound && grid_v[idx].x > 0)
                    grid_v[idx].x = 0;

                if (j < bound && grid_v[idx].y < 0)
                    grid_v[idx].y = 0;
                if (j > Scene::numGrid - bound && grid_v[idx].y > 0)
                    grid_v[idx].y = 0;

                if (k < bound && grid_v[idx].z < 0)
                    grid_v[idx].z = 0;
                if (k > Scene::numGrid - bound && grid_v[idx].z > 0)
                    grid_v[idx].z = 0;
            }
        }
        grid.sync();

        // g2p
        for (int i = 0; i < pgRepeatTime; i++) {
            size_t idx = threadId * pgRepeatTime + i;
            if (idx > Scene::numParticles) continue;
            auto base = (ivec3) (particles[idx].position * Scene::inv_dx - 0.5f);
            auto fx = (particles[idx].position * Scene::inv_dx) - (vec3) base;
            // quadratic B-spline weights
            vec3 w[] = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                        0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                        0.5f * (fx - 0.5f) * (fx - 0.5f)};
            auto new_v = vec3(0);
            auto new_C = mat3(0);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        auto offset = ivec3(i, j, k);
                        auto index = base + offset;
                        if (!(index[0] < Scene::numGrid && index[1] < Scene::numGrid &&
                              index[2] < Scene::numGrid))
                            continue;
                        auto weight = w[i][0] * w[j][1] * w[k][2];
                        auto dpos = ((vec3) offset - fx);
                        auto g_v = grid_v[index[0] * Scene::numGrid * Scene::numGrid + index[1] * Scene::numGrid + index[2]];
                        new_v += weight * g_v;
                        new_C += 4 * weight * outerProduct(g_v, dpos) * Scene::inv_dx;
                    }
                }
            }

            particles[idx].velocity = new_v;
            particles[idx].position += particles[idx].velocity * Scene::dt;
            particles[idx].J *= 1 + Scene::dt * (new_C[0][0] + new_C[1][1] + new_C[2][2]); // trace -- scale of volume
            particles[idx].C = new_C;
        }
        grid.sync();
    }
}

void Scene::update() {
    gpuUpdate();
}

void Scene::gpuUpdate() {
    cudaMemcpy(particles_gpu, &particles[0], particles_size, cudaMemcpyHostToDevice);
//    gpuCompute<<<numBlock, numThread>>>(particles_gpu, grid_v_gpu, grid_m_gpu);
    dim3 dimBlock(numThread, 1, 1);
    dim3 dimGrid(numBlock, 1, 1);
    void *kernelArgs[] = {
            (void *) &particles_gpu, (void *) &grid_v_gpu, (void *) &grid_m_gpu,
    };
    cudaError_t code = cudaLaunchCooperativeKernel((void *) gpuCompute, dimGrid, dimBlock, kernelArgs);
    if (code != cudaSuccess) {
        fprintf(stderr, "GPU assert: %s %s %d\n", cudaGetErrorString(code), __FILE__, __LINE__);
        exit(code);
    }
    cudaDeviceSynchronize();
    cudaMemcpy(&particles[0], particles_gpu, particles_size, cudaMemcpyDeviceToHost);
}

void Scene::p2g() {
    // p2g
    for (auto &particle : particles) {
        auto base = (ivec3) (particle.position * inv_dx - 0.5f);
        auto fx = (particle.position * inv_dx) - (vec3) base;
        // quadratic B-spline weights
        vector<vec3> w = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                          0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                          0.5f * (fx - 0.5f) * (fx - 0.5f)};
        Real stress = -dt * p_vol * (particle.J - 1) * 4 * inv_dx * inv_dx * E;
        mat3 affine = mat3(vec3(stress, 0, 0), vec3(0, stress, 0), vec3(0, 0, stress)) + p_mass * particle.C;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    auto offset = ivec3(i, j, k);
                    auto dpos = ((vec3) offset - fx) * dx;
                    auto weight = w[i][0] * w[j][1] * w[k][2];
                    auto index = base + offset;
                    assert (index[0] < numGrid && index[1] < numGrid && index[2] < numGrid);
                    grid_v[index[0]][index[1]][index[2]] += weight * (p_mass * particle.velocity + affine * dpos);
                    grid_m[index[0]][index[1]][index[2]] += weight * p_mass;
                }
            }
        }
    }
}

void Scene::g2p() {
    // g2p
    for (auto &particle : particles) {
        auto base = (ivec3) (particle.position * inv_dx - 0.5f);
        auto fx = (particle.position * inv_dx) - (vec3) base;
        // quadratic B-spline weights
        vector<vec3> w = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                          0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                          0.5f * (fx - 0.5f) * (fx - 0.5f)};

        auto new_v = vec3(0);
        auto new_C = mat3(0);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    auto offset = ivec3(i, j, k);
                    auto index = base + offset;
                    assert (index[0] < numGrid && index[1] < numGrid && index[2] < numGrid);
                    auto weight = w[i][0] * w[j][1] * w[k][2];
                    auto dpos = ((vec3) offset - fx) * dx;
                    auto g_v = grid_v[index[0]][index[1]][index[2]];
                    new_v += weight * g_v;
                    new_C += 4 * weight * outerProduct(g_v, dpos) * inv_dx;
                }
            }
        }

        particle.velocity = new_v;
        particle.position += particle.velocity * dt;
        particle.J *= 1 + dt * (new_C[0][0] + new_C[1][1] + new_C[2][2]); // trace -- scale of volume
        particle.C = new_C;
    }
}

void Scene::gridCompute() {
    for (int i = 0; i < numGrid; i++) {
        for (int j = 0; j < numGrid; j++) {
            for (int k = 0; k < numGrid; k++) {
                if (grid_m[i][j][k] > 0) {
                    auto inv_m = 1.0f / grid_m[i][j][k];
                    grid_v[i][j][k] = inv_m * grid_v[i][j][k];
                    grid_v[i][j][k][1] -= dt * 9.8;
                    auto bound = 3;

                    if (i < bound && grid_v[i][j][k].x < 0)
                        grid_v[i][j][k].x = 0;
                    if (i > numGrid - bound && grid_v[i][j][k].x > 0)
                        grid_v[i][j][k].x = 0;

                    if (j < bound && grid_v[i][j][k].y < 0)
                        grid_v[i][j][k].y = 0;
                    if (j > numGrid - bound && grid_v[i][j][k].y > 0)
                        grid_v[i][j][k].y = 0;

                    if (k < bound && grid_v[i][j][k].z < 0)
                        grid_v[i][j][k].z = 0;
                    if (k > numGrid - bound && grid_v[i][j][k].z > 0)
                        grid_v[i][j][k].z = 0;
                }
            }
        }
    }
}

void Scene::subStep() {
    p2g();
    gridCompute();
    g2p();
}

void Scene::gpuInit() {
    cudaMalloc((void **) &particles_gpu, particles_size);
    cudaMalloc((void **) &grid_v_gpu, grid_v_size);
    cudaMalloc((void **) &grid_m_gpu, grid_m_size);
}

void Scene::gpuFree() {
    cudaFree(particles_gpu);
    cudaFree(grid_v_gpu);
    cudaFree(grid_m_gpu);
}

void Scene::render() {
    shader.use();
    shader.set("projection",
               glm::perspective(glm::radians(camera.Zoom), (float) SCR_WIDTH / (float) SCR_HEIGHT, 0.1f, 100.0f));
    shader.set("view", camera.getViewMatrix());
    glBindVertexArray(this->VAO);
    for (auto &particle : this->particles) {
        mat4 model = mat4(1.0f); // make sure to initialize matrix to identity matrix first
        model = glm::translate(model, particle.position);
        shader.set("model", model);
        shader.set("color", vec4(0, 0.5, 1, 1));
        glDrawArrays(GL_TRIANGLES, 0, 36);
    }
    glBindVertexArray(0);
}

Scene::Scene() : VAO(0),
                 grid_v(vector<vector<vector<vec3> > >
                                (numGrid, vector<vector<vec3>>(
                                        numGrid, vector<vec3>(
                                                numGrid, vec3())))),
                 grid_m(vector<vector<vector<Real> > >
                                (numGrid, vector<vector<Real>>(
                                        numGrid, vector<Real>(
                                                numGrid, 0)))),
                 particles(vector<Particle>(numParticles)) {
    gpuInit();

    GLuint VBO;
    GLfloat particleQuad[] = {
            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            -0.5f, 0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,

            -0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f, 0.5f,

            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,

            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            -0.5f, -0.5f, 0.5f,
            -0.5f, -0.5f, -0.5f,

            -0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
    };
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(particleQuad), particleQuad, GL_STATIC_DRAW);

    // Set mesh attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *) nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    camera = Camera(glm::vec3(0.25f, 0.9f, 1.7f), glm::vec3(0.0f, 1.0f, 0.0f), -75, -30);
}

Scene::~Scene() {
    glDeleteVertexArrays(1, &VAO);
    gpuFree();
}

void Scene::processInput(GLFWwindow *window, double deltaTime) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.processKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.processKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.processKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.processKeyboard(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS)
        update();
}

void Scene::loadShader(const GLchar *vertexPath, const GLchar *fragmentPath) {
    shader.load(vertexPath, fragmentPath);
}

void Scene::processMouseMovement(float xOffset, float yOffset) {
    camera.processMouseMovement(xOffset, yOffset);
}

void Scene::processMouseScroll(float yOffset) {
    camera.processMouseScroll(yOffset);
}
