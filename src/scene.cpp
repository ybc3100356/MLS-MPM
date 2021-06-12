//
// Created by ybc on 2021/6/11.
//

#include "scene.h"
#include <cassert>

void Scene::update() {
    for (int i = 0; i < 10; i++) {
        for (auto &m : grid_m)
            std::fill(m.begin(), m.end(), 0);
        for (auto &v : grid_v)
            std::fill(v.begin(), v.end(), vec2(0, 0));
        subStep();
    }
}

void Scene::subStep() {
    // p2g
    for (auto &particle : particles) {
        auto base = (ivec2) (particle.position * inv_dx - 0.5f);
        auto fx = (particle.position * inv_dx) - (vec2) base;
        // quadratic B-spline weights
        vector<vec2> w = {0.5f * ((1.5f - fx) * (1.5f - fx)),
                          0.75f - ((fx - 1.0f) * (fx - 1.0f)),
                          0.5f * (fx - 0.5f) * (fx - 0.5f)};
        Real stress = -dt * p_vol * (particle.J - 1) * 4 * inv_dx * inv_dx * E;
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
        particle.J *= 1 + dt * (new_C[0][0] + new_C[1][1]); // trace -- scale of volume
        particle.C = new_C;
    }
}

void Scene::render() {
    shader.use();
    for (auto &particle : this->particles) {
        shader.set("offset", particle.position);
        shader.set("color", vec4(0, 0.5, 1, 1));
        glBindVertexArray(this->VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    }
}

Scene::Scene(const Shader &shader) : shader(const_cast<Shader &>(shader)), VAO(0),
                                     grid_v(vector<vector<vec2>>(numGrid, vector<vec2>(numGrid, vec2(0, 0)))),
                                     grid_m(vector<vector<Real>>(numGrid, vector<Real>(numGrid, 0))),
                                     particles(vector<Particle>(numParticles)) {
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
    glDeleteVertexArrays(1, &VAO);
}

