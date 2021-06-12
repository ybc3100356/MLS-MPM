#version 330 core
layout (location = 0) in vec2 vertex;// <vec2 position>

out vec4 ParticleColor;

uniform vec2 offset;
uniform vec4 color;

void main()
{
    float scale = 0.01f;
    ParticleColor = color;
    gl_Position = vec4(((vertex.xy) * scale) + offset - vec2(1, 1), 0.0, 1.0);
}