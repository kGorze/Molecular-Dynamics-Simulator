#pragma once
#ifndef RENDER_H
#define RENDER_H


#ifndef GLAD_H
#define GLAD_H
#include <glad/glad.h>
#endif



#ifndef GLFW_H
#define GLFW_H
#include <GLFW/glfw3.h>
#endif

#ifndef GLM_H
#define GLM_H
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#endif

#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void buildCircle(float radius, int vCount, std::vector<glm::vec3>* vertices, std::vector<unsigned int>* indices);
GLFWwindow* initializeScreen(unsigned int SCR_WIDTH, unsigned int SCR_HEIGHT);
void renderAtoms(GLFWwindow* window, std::vector<std::vector<std::vector<double>>>& dataCoords);
void terminateScreen();





#endif // !RENDER_H