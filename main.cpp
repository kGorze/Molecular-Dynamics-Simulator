#include "Headers/director.h"
#include "Headers/builder.h"
#include "Headers/command.h"
#include "Headers/concretecommands.h"
#include "Headers/menu.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include <thread>
#include <atomic>

int main() {
    Director director;
    auto builder2D = std::make_shared<Simulation2DBuilder>();
    director.set_builder(builder2D);

    //CODE FOR PATH
    std::string configFilePath = R"(C:\Users\konrad_guest\CLionProjects\MDS\Parameters\config.ini)";
    auto simulation2D = std::dynamic_pointer_cast<Simulation2D>(builder2D->getSimulation());
    if (simulation2D) {
        director.initializeSimulationWithConfig(configFilePath);
        director.buildMinimalSimulation();
        director.getBuilder()->getSimulation()->printConfig();
        director.getBuilder()->getSimulation()->run(1);

        //getting the coordinates
        director.getBuilder()->openCoordinatesFile("Resources/results.csv", 1);
        director.getBuilder()->writeCoordinatesToFile(simulation2D->getDataCoordinates());
        director.getBuilder()->closeCoordinatesFile();

        director.getBuilder()->getSimulation()->setSummaryIteration();
        std::this_thread::sleep_for(std::chrono::seconds(2));
        director.getBuilder()->getSimulation()->printSummary();
        //simulation2D->printAtomsSize();
        //simulation2D->printAtomsSize();
        //simulation2D->printAtoms();
    }

    return 0;


}