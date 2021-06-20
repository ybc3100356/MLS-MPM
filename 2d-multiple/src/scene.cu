//
// Created by ybc on 2021/6/11.
//

#include "scene.h"
#include <cassert>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cooperative_groups.h>

void Scene::update() {
    gpuUpdate();
}

const int blockNum = 32;
const int threadNum = 128;
constexpr Real e = 2.7182818284590452f;

__global__ void gpuCompute(Particle *particles, vec2 *grid_v, Real *grid_m) {
    cooperative_groups::grid_group grid = cooperative_groups::this_grid();
    int threadId = int(grid.thread_rank());
    int totalThreadNum = int(grid.size());

    for (int step = 0; step < Scene::steps; step++) {
        // memset
        if (threadId == 0) {
            memset(grid_v, 0, Scene::grid_v_size);
            memset(grid_m, 0, Scene::grid_m_size);
        }
        grid.sync();

        // p2g
        assert(Scene::numParticles % totalThreadNum == 0);
        int pgRepeatTime = int(Scene::numParticles / totalThreadNum);
        for (int i = 0; i < pgRepeatTime; i++) {
            size_t idx = threadId * pgRepeatTime + i;
            auto &p = particles[idx];
            auto base = (ivec2) (p.position * Scene::inv_dx - 0.5f);
            auto fx = (p.position * Scene::inv_dx) - (vec2) base;
            // quadratic B-spline weights
            vec2 w[] = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                        0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                        0.5f * (fx - 0.5f) * (fx - 0.5f)};
            p.F = (mat2(1) + Scene::dt * p.C) * p.F;
            Real h = max(0.1f, min(5.0f, pow(e, 10 * (1.0f - p.Jp)))); // hardness
            if (p.material == Particle::Jelly) {
                h = 4;
            }
            Real mu = Scene::mu_0 * h;
            Real la = Scene::lambda_0 * h;
            if (p.material == Particle::Liquid) {
                mu = 0.0f;
            }

            //svd
            mat2 U, sig, V;
            mat2 R, S;
            Real x = p.F[0][0] + p.F[1][1], y = p.F[0][1] - p.F[1][0];
            Real scale = (1.0f / sqrt(x * x + y * y));
            Real cosx = x * scale;
            Real sinx = y * scale;
            R = mat2(cosx, sinx, -sinx, cosx);
            S = transpose(R) * p.F;
            Real c = 0, s = 0;
            Real s1 = 0, s2 = 0;
            if (fabs(S[1][0]) < 1e-5f) {
                c = 1;
                s = 0;
                s1 = S[0][0];
                s2 = S[1][1];
            } else {
                Real tao = 0.5f * (S[0][0] - S[1][1]);
                Real w = sqrt(tao * tao + S[1][0] * S[1][0]);
                Real t = 0;
                if (tao > 0) {
                    t = S[1][0] / (tao + w);
                } else {
                    t = S[1][0] / (tao - w);
                }
                c = 1.0f / sqrt(t * t + 1);
                s = -t * c;
                s1 = c * c * S[0][0] - 2 * c * s * S[1][0] + s * s * S[1][1];
                s2 = s * s * S[0][0] + 2 * c * s * S[1][0] + c * c * S[1][1];
            }
            if (s1 < s2) {
                Real tmp = s1;
                s1 = s2;
                s2 = tmp;
                V = mat2(-s, -c, c, -s);
            } else {
                V = mat2(c, -s, s, c);
            }
            U = R * V;
            sig = mat2(s1, 0, 0, s2);

            Real J = 1.0f;
            for (int d = 0; d < 2; d++) {
                auto newSig = sig[d][d];
                if (p.material == Particle::Snow) {
                    newSig = min(max(sig[d][d], 1 - 2.5e-2f), 1 + 1.5e-3f);
                }
                p.Jp *= sig[d][d] / newSig;
                sig[d][d] = newSig;
                J *= newSig;
            }
            if (p.material == Particle::Liquid) {
                p.F = mat2(1) * sqrt(J);
            } else if (p.material == Particle::Snow) {
                p.F = U * sig * transpose(V);
            }
            auto stress = 2 * mu * (p.F - U * transpose(V)) * transpose(p.F) + mat2(1) * la * J * (J - 1);
            stress = stress * (-Scene::dt * Scene::p_vol * 4 * Scene::inv_dx * Scene::inv_dx);
            auto affine = stress + Scene::p_mass * p.C;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    auto offset = ivec2(i, j);
                    auto dpos = ((vec2) offset - fx) * Scene::dx;
                    auto weight = w[i][0] * w[j][1];
                    auto index = base + offset;
                    if (!(index[0] < Scene::numGrid && index[1] < Scene::numGrid)) continue;
                    auto dv = weight * (Scene::p_mass * p.velocity + affine * dpos);
                    auto target_idx = index[0] * Scene::numGrid + index[1];
                    atomicAdd(&(grid_v[target_idx][0]), dv[0]);
                    atomicAdd(&(grid_v[target_idx][1]), dv[1]);
                    atomicAdd(&(grid_m[target_idx]), weight * Scene::p_mass);
                }
            }
        }
        grid.sync();

        // grid
        assert((Scene::numGrid * Scene::numGrid) % totalThreadNum == 0);
        int gridRepeatTime = int((Scene::numGrid * Scene::numGrid) / totalThreadNum);
        for (int i = 0; i < gridRepeatTime; i++) {
            size_t idx = threadId * gridRepeatTime + i;
            if (grid_m[idx] > 0) {
                auto inv_m = 1.0f / grid_m[idx];
                grid_v[idx] = inv_m * grid_v[idx];
                grid_v[idx][1] -= Scene::dt * 30;
                auto bound = 3;
                size_t i = idx / Scene::numGrid;
                size_t j = idx % Scene::numGrid;
                if (i < bound && grid_v[idx][0] < 0)
                    grid_v[idx][0] = 0;
                if (i > Scene::numGrid - bound && grid_v[idx][0] > 0)
                    grid_v[idx][0] = 0;
                if (j < bound && grid_v[idx][1] < 0)
                    grid_v[idx][1] = 0;
                if (j > Scene::numGrid - bound && grid_v[idx][1] > 0)
                    grid_v[idx][1] = 0;
            }
        }
        grid.sync();

        // g2p
        for (int i = 0; i < pgRepeatTime; i++) {
            size_t idx = threadId * pgRepeatTime + i;
            auto base = (ivec2) (particles[idx].position * Scene::inv_dx - 0.5f);
            auto fx = (particles[idx].position * Scene::inv_dx) - (vec2) base;
            // quadratic B-spline weights
            vec2 w[] = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                        0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                        0.5f * (fx - 0.5f) * (fx - 0.5f)};
            auto new_v = vec2(0, 0);
            auto new_C = mat2(0);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    auto offset = ivec2(i, j);
                    auto index = base + offset;
                    if (!(index[0] < Scene::numGrid && index[1] < Scene::numGrid)) continue;
                    auto weight = w[i][0] * w[j][1];
                    auto dpos = ((vec2) offset - fx) * Scene::dx;
                    auto g_v = grid_v[index[0] * Scene::numGrid + index[1]];
                    new_v += weight * g_v;
                    new_C += 4 * weight * outerProduct(g_v, dpos) * Scene::inv_dx;
                }
            }

            particles[idx].velocity = new_v;
            particles[idx].C = new_C;
            particles[idx].position += particles[idx].velocity * Scene::dt; // boundary
        }
        grid.sync();
    }
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

void Scene::gpuUpdate() {
    cudaMemcpy(particles_gpu, &particles[0], particles_size, cudaMemcpyHostToDevice);
//    gpuCompute<<<blockNum, threadNum>>>(particles_gpu, grid_v_gpu, grid_m_gpu);
//    cudaError_t code = cudaPeekAtLastError();
    dim3 dimBlock(threadNum, 1, 1);
    dim3 dimGrid(blockNum, 1, 1);
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

void Scene::subStep() {
    p2g();
    gridCompute();
    g2p();
}

void Scene::render() {
    shader.use();
    for (auto &particle : this->particles) {
        shader.set("offset", particle.position);
        shader.set("color", particle.color);
        glBindVertexArray(this->VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    }
}

Scene::Scene(const Shader &shader) : shader(const_cast<Shader &>(shader)), VAO(0),
                                     grid_v(vector<vector<vec2>>(numGrid, vector<vec2>(numGrid, vec2(0, 0)))),
                                     grid_m(vector<vector<Real>>(numGrid, vector<Real>(numGrid, 0))) {
    for (int i = 0; i < numParticlesPerObject; i++)
        particles.emplace_back(vec2(0.35, 0.45), vec4(237 / 255., 85 / 255., 59 / 255., 1),
                               Particle::Jelly);
    for (int i = 0; i < numParticlesPerObject; i++)
        particles.emplace_back(vec2(0.45, 0.65), vec4(242 / 255., 177 / 255., 52 / 255., 1),
                               Particle::Liquid);
    for (int i = 0; i < numParticlesPerObject; i++)
        particles.emplace_back(vec2(0.55, 0.85), vec4(6 / 255., 133 / 255., 135 / 255., 1),
                               Particle::Snow);

    gpuInit();

    GLuint VBO;
    GLfloat particleQuad[] = {
            0.0f, 1.0f,
            1.0f, 0.0f,
            0.0f, 0.0f,

            0.0f, 1.0f,
            1.0f, 1.0f,
            1.0f, 0.0f,
    };
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(particleQuad), particleQuad, GL_STATIC_DRAW);

    // Set mesh attributes
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid *) nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

Scene::~Scene() {
//    glDeleteVertexArrays(1, &VAO);
    gpuFree();
}

void Scene::p2g() {
    // p2g
    for (auto &particle : particles) {
        auto base = (ivec2) (particle.position * inv_dx - 0.5f);
        auto fx = (particle.position * inv_dx) - (vec2) base;
        // quadratic B-spline weights
        vector<vec2> w = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                          0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                          0.5f * (fx - 0.5f) * (fx - 0.5f)};
        Real stress = -dt * p_vol * (particle.Jp - 1) * 4 * inv_dx * inv_dx * E;
        auto affine = mat2(vec2(stress, 0), vec2(0, stress)) + p_mass * particle.C;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto offset = ivec2(i, j);
                auto dpos = ((vec2) offset - fx) * dx;
                auto weight = w[i][0] * w[j][1];
                auto index = base + offset;
                assert (index[0] < numGrid && index[1] < numGrid);
                grid_v[index[0]][index[1]] += weight * (p_mass * particle.velocity + affine * dpos);
                grid_m[index[0]][index[1]] += weight * p_mass;
            }
        }
    }
}

void Scene::g2p() {
    // g2p
    for (auto &particle : particles) {
        auto base = (ivec2) (particle.position * inv_dx - 0.5f);
        auto fx = (particle.position * inv_dx) - (vec2) base;
        // quadratic B-spline weights
        vector<vec2> w = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                          0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                          0.5f * (fx - 0.5f) * (fx - 0.5f)};

        auto new_v = vec2(0, 0);
        auto new_C = mat2();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto offset = ivec2(i, j);
                auto index = base + offset;
                assert (index[0] < numGrid && index[1] < numGrid);
                auto weight = w[i][0] * w[j][1];
                auto dpos = ((vec2) offset - fx) * dx;
                auto g_v = grid_v[index[0]][index[1]];
                new_v += weight * g_v;
                new_C += 4 * weight * outerProduct(g_v, dpos) * inv_dx;
            }
        }

        particle.velocity = new_v;
        particle.position += particle.velocity * dt; // boundary
        particle.Jp *= 1 + dt * (new_C[0][0] + new_C[1][1]); // trace -- scale of volume
        particle.C = new_C;
    }
}

void Scene::gridCompute() {
    for (int i = 0; i < grid_m.size(); i++) {
        for (int j = 0; j < grid_m[0].size(); j++) {
            if (grid_m[i][j] > 0) {
                auto inv_m = 1.0f / grid_m[i][j];
                grid_v[i][j] = inv_m * grid_v[i][j];
                grid_v[i][j][1] -= dt * 9.8;
                auto bound = 3;
                if (i < bound && grid_v[i][j][0] < 0)
                    grid_v[i][j][0] = 0;
                if (i > numGrid - bound && grid_v[i][j][0] > 0)
                    grid_v[i][j][0] = 0;
                if (j < bound && grid_v[i][j][1] < 0)
                    grid_v[i][j][1] = 0;
                if (j > numGrid - bound && grid_v[i][j][1] > 0)
                    grid_v[i][j][1] = 0;
            }
        }
    }
}
