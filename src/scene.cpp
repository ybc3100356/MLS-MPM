//
// Created by ybc on 2021/6/11.
//

#include "scene.h"
#include <cassert>

void Scene::update() {
    for (int i = 0; i < 10; i++) {
        for (auto &m1 : grid_m)
            for (auto &m2 : m1)
                std::fill(m2.begin(), m2.end(), 0);
        for (auto &v1 : grid_v)
            for (auto &v2 : v1)
                std::fill(v2.begin(), v2.end(), vec3(0));
        subStep();
    }
}

void Scene::subStep() {
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
//        shader.set("offset", glm::vec2(particle.position[0], particle.position[1]));
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

    camera = Camera(glm::vec3(0.0f, 0.0f, 3.0f));
}

Scene::~Scene() {
    glDeleteVertexArrays(1, &VAO);
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

void Scene::processMouseMovement(float xoffset, float yoffset) {
    camera.processMouseMovement(xoffset, yoffset);
}

void Scene::processMouseScroll(float yoffset) {
    camera.processMouseScroll(yoffset);
}
