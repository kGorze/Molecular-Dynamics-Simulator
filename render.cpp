#pragma once
#include "shader.h"
#include "render.h"

void framebuffer_size_callback(GLFWwindow* window, int widht, int height) {
    glViewport(0, 0, widht, height);
}

void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }
}

GLFWwindow* initializeScreen(unsigned int SCR_WIDTH, unsigned int SCR_HEIGHT) {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Test", NULL, NULL);
    if (window == NULL) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Load OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
    }

    // Set up viewport and shaders
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    return window;
};


void buildCircle(float radius, int vCount, std::vector<glm::vec3>* vertices, std::vector<unsigned int>* indices)
{
    float angle = 360.0f / vCount;

    int triangleCount = vCount - 2;

    std::vector<glm::vec3> temp;
    // positions
    for (int i = 0; i < vCount; i++)
    {
        float currentAngle = angle * i;
        float x = radius * cos(glm::radians(currentAngle));
        float y = radius * sin(glm::radians(currentAngle));
        float z = 0.0f;

        vertices->push_back(glm::vec3(x, y, z));
    }

    // push indexes of each triangle points
    for (int i = 0; i < triangleCount; i++)
    {
        indices->push_back(0);
        indices->push_back(i + 1);
        indices->push_back(i + 2);
    }
};


void renderAtoms(GLFWwindow* window, std::vector<std::vector<std::vector<double>>>& dataCoords) {
    // Determine the maximum and minimum values in both x and y dimensions
    double maxX = -std::numeric_limits<double>::infinity();
    double minX = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();

    for (const auto& iteration : dataCoords) {
        for (const auto& coord : iteration) {
            double x = coord[0];
            double y = coord[1];
            maxX = std::max(maxX, x);
            minX = std::min(minX, x);
            maxY = std::max(maxY, y);
            minY = std::min(minY, y);
        }
    }

    // Calculate the ranges
    double rangeX = maxX - minX;
    double rangeY = maxY - minY;

    // Initialize GLFW and create window
    Shader ourShader("vertex.vs", "fragment1.fs");

    // Generate circle geometry
    std::vector<glm::vec3> vertices;
    std::vector<unsigned int> indices;
    buildCircle(0.01f, 32, &vertices, &indices);

    // Create and bind vertex array object and buffers
    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), vertices.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Main rendering loop
    int framesCount = -1;
    auto startTime = std::chrono::high_resolution_clock::now();
    double frameDuration = 1.0 / 100; // Targeting 100 frames per second

    std::cout << "\n-----------------------------------------------------------------------\n";
    std::cout << "Data visualization\n";
    std::cout << "-----------------------------------------------------------------------\n\n";

    while (!glfwWindowShouldClose(window)) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        double deltaTime = std::chrono::duration<double>(currentTime - startTime).count();

        // Check if it's time to render the next frame
        if (deltaTime >= frameDuration) {
            framesCount++;
            startTime = currentTime;

            processInput(window);

            // Clear the screen
            glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);

            // Render circles for each line in dataCoords
            if (framesCount % 100 == 0) {
                std::cout << "Frame: " << framesCount << std::endl;
            }

            const auto& iteration = dataCoords[(framesCount)];
            for (const auto& coord : iteration) {
                // Normalize coordinates to the range [-1, 1]
                double normalizedX = 2.0 * ((coord[0] - minX) / rangeX) - 1.0;
                double normalizedY = 2.0 * ((coord[1] - minY) / rangeY) - 1.0;

                // Create a translation matrix based on the normalized coordinates
                glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(normalizedX, normalizedY, 0.0f));

                // Apply the transformation matrix to the shader
                ourShader.use();
                ourShader.setMat4("model", model);

                // Render the circle
                glBindVertexArray(VAO);
                glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
            }

            if (framesCount == (dataCoords.size()-1)) {
                std::cout << "Frame: " << framesCount+1 << std::endl;
                break;
            }
            // Swap buffers and poll events
            glfwSwapBuffers(window);
            glfwPollEvents();
        }

    }

    // Clean up resources
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}



void terminateScreen() {
	glfwTerminate();
};
