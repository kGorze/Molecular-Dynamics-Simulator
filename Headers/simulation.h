//
// Created by konrad_guest on 16/05/2024.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <utility>
#include <iostream>
#include "abstract_factory.h"

class Simulation {
public:
    // Lifecycle methods
    virtual         ~Simulation()                                                   = default;
    // Simulation control methods
    virtual void    run(unsigned int option)                                        = 0;
    virtual double  singleSimulationStep()                                          = 0;
    // Configuration methods
    virtual void    printConfig()                                                   = 0;
    virtual void    setConfig(std::unordered_map<std::string, std::string> config)  = 0;
    // Physics methods
    virtual void    leapfrogStep(unsigned int part)                                 = 0;
    virtual void    applyBoundaryConditions()                                       = 0;
    virtual void    computeForces()                                                 = 0;

    // Property evaluation methods
    virtual void    evaluateProperties()                                            = 0;
    virtual void    evaluateVelocityDistribution()                                  = 0;
    virtual void    accumulateProperties(unsigned int code)                         = 0;
    virtual void    setIterationProperties()                                        = 0;

    // Output methods
    virtual void    printSummary() const                                            = 0;
    virtual void    printVelocityDestribution() const                               = 0;

};



#endif //SIMULATION_H
