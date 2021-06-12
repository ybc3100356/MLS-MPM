//
// Created by ybc on 2021/6/12.
//

#ifndef GL_SHADER_H
#define GL_SHADER_H

#include <glad/glad.h>
#include <glm/glm.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader {
public:
    unsigned int programID;

    Shader(const GLchar *vertexPath, const GLchar *fragmentPath);

    ~Shader() { if (programID) glDeleteProgram(programID); }

    void use();

    void set(const std::string &name, const glm::vec2 &value) const;

    void set(const std::string &name, const glm::vec4 &value) const;

    void set(const std::string &name, bool value) const;

    void set(const std::string &name, int value) const;

    void set(const std::string &name, float value) const;

private:
    static void checkCompileErrors(unsigned int shader, std::string type);
};

#endif //GL_SHADER_H
