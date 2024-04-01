#pragma once

#include "simulation.h"
#include "render.h"

void Showcase(std::vector<std::vector<std::vector<double>>> mdataCoordinates);

int main() {    
    Simulation simtemp;
    simtemp.setFiles("data.in", "coordinates.txt", "summary.txt", "histo.txt");
    simtemp.runSimulation(1);

    Showcase(simtemp.getCoordinatesData());

    return 0;
    };

void Showcase(std::vector<std::vector<std::vector<double>>> mdataCoordinates) {
    GLFWwindow* window = initializeScreen(800, 800);
    renderAtoms(window, mdataCoordinates);
    terminateScreen();
}