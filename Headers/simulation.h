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
    virtual void run(unsigned int option) = 0;
    virtual ~Simulation() = default;
    virtual void printConfig() = 0;
    virtual void leapfrogStep(unsigned int part) = 0;
    virtual double singleSimulationStep() = 0;
    virtual void applyBoundaryConditions() = 0;
    virtual void computeForcesPairs() = 0;
    virtual void evaluateProperties() = 0;
    virtual void evaluateVelocityDistribution() = 0;
    virtual void accumulateProperties(unsigned int code) = 0;
    virtual void printSummary() const = 0;;
    virtual void printVelocityDestribution() const = 0;

    virtual void setIterationProperties() = 0;
};



#endif //SIMULATION_H
