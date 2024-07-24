#include "Headers/director.h"
#include "Headers/builder.h"
#include "Headers/simulation2d.h"
#include "Headers/command.h"
#include "Headers/concretecommands.h"
#include "Headers/menu.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include <thread>
#include <atomic>

//PLOTTING
#include <matplot/matplot.h>
#include <thread>
#include <chrono>


int main(int argc, char** argv){

    // if(argc != 2) {
    //     std::cout<<"Usage: "<<argv[0]<<" <config_file_path>\n";
    //     return 1;
    // }

     Director director;
     auto builder2D = std::make_shared<Simulation2DBuilder>();
     director.set_builder(builder2D);

     //CODE FOR PATH
     std::string configFilePath = R"(C:\Users\konrad_guest\source\repos\MDS\Parameters\config.ini)";
     auto simulation2D = std::dynamic_pointer_cast<Simulation2D>(builder2D->getSimulation());
     if (simulation2D) {
         director.initializeSimulationWithConfig(configFilePath);
         director.buildMinimalSimulation();
         director.getBuilder()->getSimulation()->printConfig();
         director.getBuilder()->getSimulation()->run(1);


         builder2D->saveDataToFile("Resources/results.csv", simulation2D->getDataCoordinates(), 1);
         builder2D->saveDataToFile("Resources/properties.csv", simulation2D->getDataProperties(), 1);
         builder2D->saveDataToFile("Resources/velocity_distribution.csv", simulation2D->getDataHistogramVelocities(), 1);


         //getting the velocity distribution

         // auto dataHistogramVelocities = simulation2D->getDataHistogramVelocities();
         // simulation2D->plotDataHistogramVelocities(dataHistogramVelocities);




         std::this_thread::sleep_for(std::chrono::seconds(2));
         director.getBuilder()->getSimulation()->printSummary();
         std::cout<<"Simulation finished\n";

    }


    return 0;

}