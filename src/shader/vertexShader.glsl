#version 330 core
layout (location = 0) in vec3 aPos;// <vec2 position>

out vec4 ParticleColor;

//uniform vec2 offset;
uniform vec4 color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    float scale = 0.003f;
    ParticleColor = color;
    gl_Position = projection * view * model * vec4((aPos * scale), 1.0f);
}