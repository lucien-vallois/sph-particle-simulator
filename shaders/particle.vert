#version 330 core

layout (location = 0) in vec3 aPos;          // Quad vertex positions
layout (location = 1) in vec3 aInstancePos;  // Instance position
layout (location = 2) in vec3 aVelocity;     // Particle velocity (for color mapping)
layout (location = 3) in vec3 aColor;        // Particle color

// Uniforms
uniform mat4 view;
uniform mat4 projection;
uniform mat4 model;
uniform float particleSize;
uniform int colorMode;

// Outputs to fragment shader
out vec3 fragColor;
out vec3 fragVelocity;
out vec3 fragPosition;

void main()
{
    // Calculate world position for this instance
    vec3 worldPos = aInstancePos;

    // Scale quad based on particle size and distance for billboarding
    vec3 cameraRight = vec3(view[0][0], view[1][0], view[2][0]);
    vec3 cameraUp = vec3(view[0][1], view[1][1], view[2][1]);

    vec3 pos = worldPos +
               cameraRight * aPos.x * particleSize +
               cameraUp * aPos.y * particleSize;

    // Transform to clip space
    gl_Position = projection * view * model * vec4(pos, 1.0);

    // Pass data to fragment shader
    fragVelocity = aVelocity;
    fragPosition = aInstancePos;
    fragColor = aColor;
}

