//
// Created by konrad_guest on 24/07/2024.
//

#ifndef SIMULATION2D_H
#define SIMULATION2D_H

#include "simulation.h"

#include "boundary_conditions.h"

#include "matplot/matplot.h"

#include "leapfrog.h"

#include "errorcodes.h"

#include "abstract_factory.h"

#include "progressbar.h"

#include <unordered_map>

#include <memory>

#include <iomanip>

#include <fstream>

#include <cmath>

#include <Eigen/Dense>

#include <filesystem>

#include <chrono>

#include <ctime>

#include <filesystem>


class Simulation2D : public Simulation
{
private:
    double      timeNow,temperature, deltaT, density, rangeVel, potentialEnergySum, viralEnergySum, velocityMagnitude, cutoffRadius, deltaVelocity, histogramSum, entropyFunction;
    int         stepCount, stepAvg, stepEquil, stepLimit, limitVel, sizeHistVel, stepVel, randSeed,numberOfDimensions, numberOfAtoms, numberOfAtomIterations, countStep, countVelocities;

    std::vector<std::shared_ptr<Atom2D>>            atoms;
    std::vector<double>                             histogramVelocities;
    std::vector<std::vector<double>>                dataHistogramVelocities;
    std::vector<std::tuple<int, double, double>>    dataCoordinates;
    std::unordered_map<std::string, std::string>    config;
    std::shared_ptr<AtomFactory>                    atomFactory;
    std::shared_ptr<Progressbar>                    progressbar;

    std::vector<std::tuple<int, double, Eigen::Vector2d, double, double, double, double, double, double>> iterationData;

    Eigen::Vector2d initUcell,velocitiesSum,region;
    Eigen::Vector3d kineticEnergy, potentialEnergy, pressure;
public:
    void    setConfig(std::unordered_map<std::string, std::string> config) {this->config = config;}
    void    printConfig() override ;
    void    leapfrogStep(unsigned int part) override;
    double  singleSimulationStep() override;
    void    applyBoundaryConditions() override;
    void    computeForcesPairs() override;
    void    evaluateProperties() override;
    void    evaluateVelocityDistribution() override;
    void    accumulateProperties(unsigned int code) override;
    void    printSummary() const override;
    void    printVelocityDestribution() const override;

    void    setIterationProperties() override;

    void setDeltaT(double deltaT) { this->deltaT = deltaT; }
    void setDensity(double density) { this->density = density; }
    void setInitUcell(int x, int y) { this->initUcell(0) = x; this->initUcell(1) = y; }
    void setStepAvg(int stepAvg) { this->stepAvg = stepAvg; }
    void setStepEquil(int stepEquil) { this->stepEquil = stepEquil; }
    void setStepLimit(int stepLimit) { this->stepLimit = stepLimit; }
    void setTemperature(double temperature) { this->temperature = temperature; }
    void setLimitVel(int limitVel) { this->limitVel = limitVel; }
    void setRangeVel(double rangeVel) { this->rangeVel = rangeVel; }
    void setSizeHistVel(int sizeHistVel) { this->sizeHistVel = sizeHistVel; }
    void setStepVel(int stepVel) { this->stepVel = stepVel; }
    void setRandSeed(int randSeed) { this->randSeed = randSeed; }
    void setNumberOfDimensions(int numberOfDimensions) { this->numberOfDimensions = numberOfDimensions;}
    void setCutoffRadius(double cutoffRadius) { this->cutoffRadius = cutoffRadius; }

    void setRegion(Eigen::Vector2d region) { this->region = region; }
    void setVelocitiesSum(Eigen::Vector2d velocitiesSum) { this->velocitiesSum = velocitiesSum; }
    void setAtomFactory(std::shared_ptr<AtomFactory> atomFactory) { this->atomFactory = atomFactory; }
    void setProgressbar(std::shared_ptr<Progressbar> progressbar) { this->progressbar = progressbar; }

    void setNumberOfAtoms(int numberOfAtoms) { this->numberOfAtoms = numberOfAtoms; }
    void setVelocityMagnitude(double velocityMagnitude) { this->velocityMagnitude = velocityMagnitude; }

    void setNumberOfAtomIterations(int numberOfAtomIterations) { this->numberOfAtomIterations = numberOfAtomIterations; }
    void setCountStep(int countStep) { this->countStep = countStep; }
    void setCountVelocities(int countVelocities) { this->countVelocities = countVelocities; }
    void setpotentialEnergySum(double potentialEnergySum) { this->potentialEnergySum = potentialEnergySum; }
    void setviralEnergySum(double viralEnergySum){ this->viralEnergySum = viralEnergySum;}
    void setDeltaVelocity(double deltaVelocity) { this->deltaVelocity = deltaVelocity; }
    void setHistogramSum(double histogramSum) { this->histogramSum = histogramSum; }
    void setEntropyFunction(double entropyFunction) { this->entropyFunction = entropyFunction; }
    void setStepCount(int stepCount) { this->stepCount = stepCount; }
    void setTimeNow(double timeNow) { this->timeNow = timeNow; }

    void setKineticEnergy(Eigen::Vector3d kineticEnergy) { this->kineticEnergy = kineticEnergy; }
    void setPotentialEnergy(Eigen::Vector3d potentialEnergy) { this->potentialEnergy = potentialEnergy; }
    void setPressure(Eigen::Vector3d pressure) { this->pressure = pressure; }
    void setHistogramVelocities(std::vector<double> histogramVelocities) { this->histogramVelocities = histogramVelocities; }
    void setDataHistogramVelocities(std::vector<std::vector<double>> dataHistogramVelocities) { this->dataHistogramVelocities = dataHistogramVelocities; }

    template<typename T>
    T get(const std::string& param) const;

    std::vector<std::shared_ptr<Atom2D>>&           getAtoms() { return atoms; }
    std::vector<double>&                            getHistogramVelocities()& {return histogramVelocities;}
    std::vector<std::tuple<int, double, double>>&   getDataCoordinates() {return dataCoordinates;}
    std::vector<std::vector<double>>&               getDataHistogramVelocities() {return dataHistogramVelocities;}
    std::shared_ptr<AtomFactory>&                   getAtomFactory() {return atomFactory;}
    std::shared_ptr<Progressbar>&                   getProgressbar() { return progressbar; }

    std::vector<std::tuple<int, double, Eigen::Vector2d, double, double, double, double, double, double>>& getDataProperties() {return iterationData;}

    static void plotDataHistogramVelocities(const std::vector<std::vector<double>>& dataHistogramVelocities);

    Simulation2D() = default;
    void run(unsigned int option) override;
};



template<>
inline double Simulation2D::get<double>(const std::string& param) const
{
    if (param == "deltaT") return deltaT;
    if (param == "density") return density;
    if (param == "rangeVel") return rangeVel;
    if (param == "temperature") return temperature;
    if (param == "velocityMagnitude") return velocityMagnitude;
    if (param == "cutOffRadius") return cutoffRadius;
    if (param == "viralEnergySum") return viralEnergySum;
    if (param == "deltaVelocity") return deltaVelocity;
    if (param == "histogramSum") return histogramSum;
    if (param == "entropyFunction") return entropyFunction;
    if (param == "potentialEnergySum") return potentialEnergySum;
    if (param == "timeNow") return timeNow;
    throw std::invalid_argument("Invalid parameter name for type double: " + param);
}

template<>
inline int Simulation2D::get<int>(const std::string& param) const
{
    if (param == "initUcellX") return static_cast<int>(initUcell(0));
    if (param == "initUcellY") return static_cast<int>(initUcell(1));
    if (param == "stepAvg") return stepAvg;
    if (param == "stepEquil") return stepEquil;
    if (param == "stepLimit") return stepLimit;
    if (param == "limitVel") return limitVel;
    if (param == "sizeHistVel") return sizeHistVel;
    if (param == "stepVel") return stepVel;
    if (param == "randSeed") return randSeed;
    if (param == "numberOfDimensions") return numberOfDimensions;
    if (param == "numberOfAtoms") return numberOfAtoms;
    if (param == "numberOfAtomIterations") return numberOfAtomIterations;
    if (param == "countStep") return countStep;
    if (param == "countVelocities") return countVelocities;
    if (param == "stepCount") return stepCount;
    throw std::invalid_argument("Invalid parameter name for type int: " + param);
}

template<>
inline Eigen::Vector2d Simulation2D::get<Eigen::Vector2d>(const std::string& param) const
{
    if (param == "region") return region;
    if (param == "initUcell") return initUcell;
    if (param == "velocitiesSum") return velocitiesSum;
    throw std::invalid_argument("Invalid parameter name for type Eigen::Vector2d: " + param);
}

template<>
inline Eigen::Vector3d Simulation2D::get<Eigen::Vector3d>(const std::string& param) const
{
    if (param == "kineticEnergy") return kineticEnergy;
    if (param == "potentialEnergy") return potentialEnergy;
    if (param == "pressure") return pressure;
    throw std::invalid_argument("Invalid parameter name for type Eigen::Vector3d: " + param);
}






#endif //SIMULATION2D_H
