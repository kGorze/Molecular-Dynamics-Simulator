//
// Created by konrad_guest on 16/05/2024.
//

#ifndef BUILDER_H
#define BUILDER_H



#include "simulation.h"
#include "abstract_factory.h"
#include "progressbar.h"
#include <unordered_map>
#include <memory>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <filesystem>


class Builder {
public:
    Builder();
    virtual ~Builder();

    virtual void setConfig(std::unordered_map<std::string, std::string> config) = 0;
    virtual void setRestOfParameters() const = 0;
    virtual void initializeVectors() = 0;
    virtual void initializeCoordinates() = 0;
    virtual void initializeVelocities() = 0;
    virtual void initializeAccelerations() = 0;
    virtual void setupStaticSimulation() = 0;

    virtual void ensureDirectoryExists(const std::string& path) = 0;


    //PARAMETERS SETTERS
    virtual void setDeltaT(double deltaT) = 0;
    virtual void setDensity(double density) = 0;
    virtual void setInitUcell(int x, int y) = 0;
    virtual void setStepAvg(int stepAvg) = 0;
    virtual void setStepEquil(int stepEquil) = 0;
    virtual void setStepLimit(int stepLimit) = 0;
    virtual void setTemperature(double temperature) = 0;
    virtual void setLimitVel(int limitVel) = 0;
    virtual void setRangeVel(double rangeVel) = 0;
    virtual void setSizeHistVel(int sizeHistVel) = 0;
    virtual void setStepVel(int stepVel) = 0;
    virtual void setRandSeed(int randSeed) = 0;
    virtual void setNumberOfDimensions(int numberOfDimensions) = 0;

    //PARAMETERS GETTERS
    virtual double getDeltaT()& = 0;
    virtual int getinitUcellX()& = 0;
    virtual int getinitUcellY()& = 0;
    virtual std::ofstream& getcoordinatesDataOutput() = 0;
    virtual std::ofstream& getPropertiesDataOutput() = 0;
    virtual std::ofstream& getVelocityDistributionDataOutput() = 0;


    virtual void initializeParameters() = 0;

    // atom simulation
    /*virtual void setAtomSize(int x_axis, int y_axis) = 0;
    virtual void initializeAtoms() = 0;
     */

    //molecule simulation
    //virtual void setMoleculeSize() const = 0;
    //virtual void initializeMolecules() const = 0;
    virtual std::shared_ptr<Simulation> getSimulation() = 0;
};

class Simulation2D : public Simulation {
private:
    //ATRIBUTES
    double deltaT, density, rangeVel;
    int stepAvg, stepEquil, stepLimit, limitVel, sizeHistVel, stepVel, randSeed;
    int numberOfDimensions;
    double cutoffRadius;
    double velocityMagnitude;
    double potentialEnergySum, viralEnergySum;
    int numberOfAtoms;
    int numberOfAtomIterations;
    int countStep, countVelocities;

    double deltaVelocity, histogramSum;
    double entropyFunction;

    int stepCount;
    double timeNow;

    std::vector<std::shared_ptr<Atom2D>> atoms;
    std::vector<double> histogramVelocities;
    std::vector<std::vector<double>> dataHistogramVelocities;
    std::vector<std::tuple<int, double, double>> dataCoordinates;

    std::unordered_map<std::string, std::string> config;
    std::shared_ptr<AtomFactory> atomFactory;
    std::shared_ptr<Progressbar> progressbar;

    std::vector<std::tuple<int, double, Eigen::Vector2d, double, double, double, double, double, double>> iterationData;


    Eigen::Vector2d region;
    Eigen::Vector2d initUcell;
    Eigen::Vector2d velocitiesSum;
    Eigen::Vector3d kineticEnergy, potentialEnergy, pressure;

    double temperature;
public:

    void setConfig(std::unordered_map<std::string, std::string> config) {
        this->config = config;
    }
    void printConfig() override ;
    void leapfrogStep(unsigned int part) override;
    double singleSimulationStep() override;
    void applyBoundaryConditions() override;
    void computeForcesPairs() override;
    void evaluateProperties() override;
    void evaluateVelocityDistribution() override;
    void accumulateProperties(unsigned int code) override;
    void printSummary() const override;
    void printVelocityDestribution() const override;

    void setIterationProperties() override;

    //PARAMETERS SETTER
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

    // PARAMETERS GETTER
    template<typename T>
    T get(const std::string& param) const;

    std::vector<std::shared_ptr<Atom2D>>& getAtoms() { return atoms; }
    std::vector<double>& getHistogramVelocities()& {return histogramVelocities;}
    std::vector<std::tuple<int, double, double>>& getDataCoordinates() {return dataCoordinates;}
    std::vector<std::vector<double>>& getDataHistogramVelocities() {return dataHistogramVelocities;}
    std::vector<std::tuple<int, double, Eigen::Vector2d, double, double, double, double, double, double>>& getDataProperties() {return iterationData;}
    std::shared_ptr<AtomFactory>& getAtomFactory() {return atomFactory;}
    std::shared_ptr<Progressbar>& getProgressbar() { return progressbar; }

    //
    static void plotDataHistogramVelocities(const std::vector<std::vector<double>>& dataHistogramVelocities);

    Simulation2D() = default;
    void run(unsigned int option) override;

    //void Reset();




    /*     void setAtomSize(int x_axis, int y_axis);
        void initializeAtoms();

        void printAtomsSize() const;
        void printAtoms() const;*/



};



template<>
inline double Simulation2D::get<double>(const std::string& param) const {
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
inline int Simulation2D::get<int>(const std::string& param) const {
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
inline Eigen::Vector2d Simulation2D::get<Eigen::Vector2d>(const std::string& param) const {
    if (param == "region") return region;
    if (param == "initUcell") return initUcell;
    if (param == "velocitiesSum") return velocitiesSum;
    throw std::invalid_argument("Invalid parameter name for type Eigen::Vector2d: " + param);
}
template<>
inline Eigen::Vector3d Simulation2D::get<Eigen::Vector3d>(const std::string& param) const {
    if (param == "kineticEnergy") return kineticEnergy;
    if (param == "potentialEnergy") return potentialEnergy;
    if (param == "pressure") return pressure;
    throw std::invalid_argument("Invalid parameter name for type Eigen::Vector3d: " + param);
}




class Simulation2DBuilder : public Builder {
private:
    std::shared_ptr<Simulation2D> simulation;
    std::ofstream coordinatesDataOutput;
    std::ofstream propertiesDataOutput;
    std::ofstream velocityDistributionDataOutput;

public:
    void setConfig(std::unordered_map<std::string, std::string> config) override {
        simulation->setConfig(config);
    }
    //PARAMETERS SETTERS
    void setDeltaT(const double deltaT) override { simulation->setDeltaT(deltaT); }
    void setDensity(const double density) override{ simulation->setDensity(density); }
    void setInitUcell(const int x,const  int y) override{ simulation->setInitUcell(x, y); }
    void setStepAvg(const int stepAvg) override{ simulation->setStepAvg(stepAvg); }
    void setStepEquil(const int stepEquil) override{ simulation->setStepEquil(stepEquil); }
    void setStepLimit(const int stepLimit) override{ simulation->setStepLimit(stepLimit); }
    void setTemperature(const double temperature) override{ simulation->setTemperature(temperature); }
    void setLimitVel(const int limitVel) override{ simulation->setLimitVel(limitVel); }
    void setRangeVel(const double rangeVel) override{ simulation->setRangeVel(rangeVel); }
    void setSizeHistVel(const int sizeHistVel) override{ simulation->setSizeHistVel(sizeHistVel); }
    void setStepVel(const int stepVel) override{ simulation->setStepVel(stepVel); }
    void setRandSeed(const int randSeed) override{ simulation->setRandSeed(randSeed); }
    void setNumberOfDimensions(int numberOfDimensions) override{ simulation->setNumberOfDimensions(numberOfDimensions);}
    void setAtomFactory(std::shared_ptr<AtomFactory> atomFactory) { simulation->setAtomFactory(atomFactory); }
    void setNumberOfAtomIterations(int numberOfAtomIterations) { simulation->setNumberOfAtomIterations(numberOfAtomIterations); }
    void setCoordinatesDataOutput(std::ofstream&& coordinatesDataOutput) {this->coordinatesDataOutput = std::move(coordinatesDataOutput);}
    void setPropertiesDataOutput(std::ofstream&& propertiesDataOutput) {this->propertiesDataOutput = std::move(propertiesDataOutput);}
    void setVelocityDistributionDataOutput(std::ofstream&& velocityDistributionDataOutput) {this->velocityDistributionDataOutput = std::move(velocityDistributionDataOutput);}

    ////PARAMETERS GETTERS
    double getDeltaT()& override { return simulation->get<double>("deltaT"); }
    int getinitUcellX()& override { return simulation->get<int>("initUcellX"); }
    int getinitUcellY()& override { return simulation->get<int>("initUcellY"); }
    std::vector<std::shared_ptr<Atom2D>>& getAtoms()  { return simulation->getAtoms();}
    std::vector<double>& getHistogramVelocities()& { return simulation->getHistogramVelocities();}
    std::vector<std::tuple<int, double, double>>& getDataCoordinates() { return simulation->getDataCoordinates();}
    std::vector<std::vector<double>>& getDataHistogramVelocities() { return simulation->getDataHistogramVelocities();}
    std::shared_ptr<AtomFactory>& getAtomFactory() {return simulation->getAtomFactory();}
    std::ofstream& getcoordinatesDataOutput() override { return coordinatesDataOutput; }
    std::ofstream& getPropertiesDataOutput() { return propertiesDataOutput; }
    std::ofstream& getVelocityDistributionDataOutput() { return velocityDistributionDataOutput; }


    Simulation2DBuilder() {
        reset();
    }

    void setRestOfParameters() const ;
    void initializeVectors() override;
    void initializeCoordinates() override;
    void initializeVelocities() override;
    void initializeAccelerations() override;
    void setupStaticSimulation() override;

    void ensureDirectoryExists(const std::string& path) override;

    template <typename T>
    void saveDataToFile(const std::string& filename, const T& data,  unsigned int mode) {
        static_assert(sizeof(T) == -1, "This function must be used with a specific specialization.");
    }

    void initializeParameters() override{
        simulation->setTemperature(10);
    }

    void reset() {
        simulation = std::make_shared<Simulation2D>();
    }
    /*     void setAtomSize(int x_axis, int y_axis) override {
            simulation->setAtomSize(x_axis, y_axis);
        }
        void initializeAtoms() override {
            simulation->initializeAtoms();
        } */
    std::shared_ptr<Simulation> getSimulation() override;
};

template <> //template for the properties
void Simulation2DBuilder::saveDataToFile(
    const std::string& filename,
    const std::vector<std::tuple<int, double, Eigen::Vector2d, double, double, double, double, double, double>>& data,
    unsigned int mode);

template <> //template for the velocities
void Simulation2DBuilder::saveDataToFile(
    const std::string& filename,
    const std::vector<std::vector<double>>& data,
    unsigned int mode);


template <> //template for the coordinates
void Simulation2DBuilder::saveDataToFile(
    const std::string& filename,
    const std::vector<std::tuple<int, double, double>>& data,
    unsigned int mode);

#endif //BUILDER_H
