#version 330 core

// Inputs from vertex shader
in vec3 fragColor;
in vec3 fragVelocity;
in vec3 fragPosition;

// Uniforms
uniform int colorMode;  // 0: velocity, 1: density, 2: pressure, 3: type

// Output
out vec4 FragColor;

// Color mapping functions
vec3 hsv2rgb(vec3 hsv) {
    vec3 rgb;
    float h = hsv.x;
    float s = hsv.y;
    float v = hsv.z;

    float c = v * s;
    float x = c * (1.0 - abs(mod(h / 60.0, 2.0) - 1.0));
    float m = v - c;

    if (h >= 0.0 && h < 60.0) {
        rgb = vec3(c, x, 0.0);
    } else if (h >= 60.0 && h < 120.0) {
        rgb = vec3(x, c, 0.0);
    } else if (h >= 120.0 && h < 180.0) {
        rgb = vec3(0.0, c, x);
    } else if (h >= 180.0 && h < 240.0) {
        rgb = vec3(0.0, x, c);
    } else if (h >= 240.0 && h < 300.0) {
        rgb = vec3(x, 0.0, c);
    } else {
        rgb = vec3(c, 0.0, x);
    }

    return rgb + vec3(m);
}

vec3 velocityToColor(vec3 velocity) {
    float speed = length(velocity);
    float maxSpeed = 2.0;  // Adjust based on simulation
    float normalizedSpeed = clamp(speed / maxSpeed, 0.0, 1.0);

    // Blue (slow) to Red (fast)
    return mix(vec3(0.0, 0.5, 1.0), vec3(1.0, 0.0, 0.0), normalizedSpeed);
}

vec3 densityToColor(float density) {
    // Assume rest density ~1000, max density ~2000
    float restDensity = 1000.0;
    float maxDensity = 2000.0;
    float normalizedDensity = clamp((density - restDensity) / (maxDensity - restDensity), 0.0, 1.0);

    // Light blue (low density) to dark blue (high density)
    return mix(vec3(0.7, 0.9, 1.0), vec3(0.0, 0.2, 0.5), normalizedDensity);
}

vec3 pressureToColor(float pressure) {
    float maxPressure = 10000.0;
    float normalizedPressure = clamp(pressure / maxPressure, 0.0, 1.0);

    // Green (low pressure) to Yellow/Red (high pressure)
    return mix(vec3(0.0, 0.8, 0.0), vec3(1.0, 0.0, 0.0), normalizedPressure);
}

void main()
{
    vec3 finalColor;

    if (colorMode == 0) {
        // Velocity-based coloring
        finalColor = velocityToColor(fragVelocity);
    } else if (colorMode == 1) {
        // Density-based coloring (would need density as input)
        // For now, use position-based coloring as placeholder
        float height = fragPosition.y;
        finalColor = densityToColor(1000.0 + height * 200.0);
    } else if (colorMode == 2) {
        // Pressure-based coloring (would need pressure as input)
        // For now, use velocity magnitude as placeholder
        finalColor = pressureToColor(length(fragVelocity) * 5000.0);
    } else {
        // Type-based coloring (use stored color)
        finalColor = fragColor;
    }

    // Add some lighting effect
    vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
    vec3 normal = normalize(vec3(0.0, 0.0, 1.0));  // Billboard normal
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 ambient = finalColor * 0.3;
    vec3 diffuse = finalColor * diff * 0.7;

    finalColor = ambient + diffuse;

    // Distance-based alpha for depth cueing
    float alpha = 1.0;
    // Could add distance-based fading here

    FragColor = vec4(finalColor, alpha);
}

