#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>

#include "scene.h"

static void processInput(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

static void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    glViewport(0, 0, width, height);
}

static void mouse_callback(GLFWwindow *window, double xpos, double ypos) {
    Scene::processMouseMovement(xpos, ypos);
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        Scene::setMouse(1);
    else
        Scene::setMouse(0);
}

#include <cstdio>

int test() {
    mat2 A, U, S, V;
    A = mat2(-3, 7, -12, -9);
    svd(A, U, S, V);
    printf("A=\n%f, %f,\n%f, %f\n\n", A[0][0], A[1][0], A[0][1], A[1][1]);
    printf("U=\n%f, %f,\n%f, %f\n\n", U[0][0], U[1][0], U[0][1], U[1][1]);
    printf("S=\n%f, %f,\n%f, %f\n\n", S[0][0], 0.0, 0.0, S[1][1]);
    printf("V=\n%f, %f,\n%f, %f\n\n", V[0][0], V[1][0], V[0][1], V[1][1]);
    exit(0);
}

int main() {
//    test();
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // glfw window creation
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "MLS-MPM", nullptr, nullptr);
    if (window == nullptr) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    glfwSetCursorPosCallback(window, mouse_callback);
//    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    double deltaTime = 0.0f;
    double lastFrame = 0.0f;
    Shader shader = Shader("../src/shader/vertexShader.glsl", "../src/shader/fragmentShader.glsl");
    Scene scene = Scene(shader);
    while (!glfwWindowShouldClose(window)) {
        double currentFrame = glfwGetTime();
        deltaTime += currentFrame - lastFrame;
        lastFrame = currentFrame;

        glfwPollEvents();
        processInput(window);
        scene.update();
        glClearColor(17 / 255., 47 / 255., 65 / 255., 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        scene.render();
        glfwSwapBuffers(window);
    }
    // Terminate GLFW, clearing any resources allocated by GLFW.
    glfwTerminate();
    std::cout << "finished" << std::endl;
    return 0;
}
