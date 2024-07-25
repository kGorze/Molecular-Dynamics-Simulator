//
// Created by konrad_guest on 16/05/2024.
//

#ifndef BUILDER_H
#define BUILDER_H


#include "simulation.h"

#include "simulation2dap.h"

#include "datasaver.h"

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




class Builder
{
public:
    Builder();

    virtual         ~Builder();

    virtual void    setConfig(std::unordered_map<std::string, std::string> config) = 0;
    virtual void    setRestOfParameters() const = 0;
    virtual void    initializeVectors() = 0;
    virtual void    initializeCoordinates() = 0;
    virtual void    initializeVelocities() = 0;
    virtual void    initializeAccelerations() = 0;
    virtual void    setupStaticSimulation() = 0;

    virtual void    ensureDirectoryExists(const std::string& path) = 0;

    virtual void    setDeltaT(double deltaT) = 0;
    virtual void    setDensity(double density) = 0;
    virtual void    setInitUcell(int x, int y) = 0;
    virtual void    setStepAvg(int stepAvg) = 0;
    virtual void    setStepEquil(int stepEquil) = 0;
    virtual void    setStepLimit(int stepLimit) = 0;
    virtual void    setTemperature(double temperature) = 0;
    virtual void    setLimitVel(int limitVel) = 0;
    virtual void    setRangeVel(double rangeVel) = 0;
    virtual void    setSizeHistVel(int sizeHistVel) = 0;
    virtual void    setStepVel(int stepVel) = 0;
    virtual void    setRandSeed(int randSeed) = 0;
    virtual void    setNumberOfDimensions(int numberOfDimensions) = 0;

    virtual double         getDeltaT()& = 0;
    virtual int            getinitUcellX()& = 0;
    virtual int            getinitUcellY()& = 0;
    virtual std::ofstream& getcoordinatesDataOutput() = 0;
    virtual std::ofstream& getPropertiesDataOutput() = 0;
    virtual std::ofstream& getVelocityDistributionDataOutput() = 0;


    virtual void                        initializeParameters() = 0;
    virtual std::shared_ptr<Simulation> getSimulation() = 0;
};



class Simulation2DBuilder : public Builder
{
private:
    std::shared_ptr<simulation2dap>   simulation;
    std::ofstream                   coordinatesDataOutput,propertiesDataOutput,velocityDistributionDataOutput;
    DataSaver2D dataSaver;
public:
    void setConfig(std::unordered_map<std::string, std::string> config) override {simulation->setConfig(config);}


    void setDeltaT(const double deltaT) override                    { simulation->setDeltaT(deltaT); }
    void setDensity(const double density) override                  { simulation->setDensity(density); }
    void setInitUcell(const int x,const  int y) override            { simulation->setInitUcell(x, y); }
    void setStepAvg(const int stepAvg) override                     { simulation->setStepAvg(stepAvg); }
    void setStepEquil(const int stepEquil) override                 { simulation->setStepEquil(stepEquil); }
    void setStepLimit(const int stepLimit) override                 { simulation->setStepLimit(stepLimit); }
    void setTemperature(const double temperature) override          { simulation->setTemperature(temperature); }
    void setLimitVel(const int limitVel) override                   { simulation->setLimitVel(limitVel); }
    void setRangeVel(const double rangeVel) override                { simulation->setRangeVel(rangeVel); }
    void setSizeHistVel(const int sizeHistVel) override             { simulation->setSizeHistVel(sizeHistVel); }
    void setStepVel(const int stepVel) override                     { simulation->setStepVel(stepVel); }
    void setRandSeed(const int randSeed) override                   { simulation->setRandSeed(randSeed); }
    void setNumberOfDimensions(int numberOfDimensions) override     { simulation->setNumberOfDimensions(numberOfDimensions);}
    void setAtomFactory(std::shared_ptr<AtomFactory> atomFactory)   { simulation->setAtomFactory(atomFactory); }
    void setNumberOfAtomIterations(int numberOfAtomIterations)      { simulation->setNumberOfAtomIterations(numberOfAtomIterations); }
    void setCoordinatesDataOutput(std::ofstream&& coordinatesDataOutput)                    {this->coordinatesDataOutput = std::move(coordinatesDataOutput);}
    void setPropertiesDataOutput(std::ofstream&& propertiesDataOutput)                      {this->propertiesDataOutput = std::move(propertiesDataOutput);}
    void setVelocityDistributionDataOutput(std::ofstream&& velocityDistributionDataOutput)  {this->velocityDistributionDataOutput = std::move(velocityDistributionDataOutput);}


    double getDeltaT()& override    { return simulation->get<double>("deltaT"); }
    int getinitUcellX()& override   { return simulation->get<int>("initUcellX"); }
    int getinitUcellY()& override   { return simulation->get<int>("initUcellY"); }
    std::vector<std::shared_ptr<Atom>>& getAtoms()                    { return simulation->getAtoms();}
    std::vector<double>& getHistogramVelocities()&                      { return simulation->getHistogramVelocities();}
    std::vector<std::tuple<int, Eigen::VectorXd>>& getDataCoordinates()  { return simulation->getDataCoordinates();}
    std::vector<std::vector<double>>& getDataHistogramVelocities()      { return simulation->getDataHistogramVelocities();}
    std::shared_ptr<AtomFactory>& getAtomFactory()                      {return simulation->getAtomFactory();}
    std::ofstream& getcoordinatesDataOutput() override                  { return coordinatesDataOutput; }
    std::ofstream& getPropertiesDataOutput()                            { return propertiesDataOutput; }
    std::ofstream& getVelocityDistributionDataOutput()                  { return velocityDistributionDataOutput; }
    DataSaver2D& getDataSaver2D()                                       { return dataSaver; }




    Simulation2DBuilder() {reset();}

    void setRestOfParameters() const ;
    void initializeVectors() override;
    void initializeCoordinates() override;
    void initializeVelocities() override;
    void initializeAccelerations() override;
    void setupStaticSimulation() override;

    void ensureDirectoryExists(const std::string& path) override;

    void initializeParameters() override{simulation->setTemperature(10);}

    void reset() {simulation = std::make_shared<simulation2dap>();}
    std::shared_ptr<Simulation> getSimulation() override;
};

#endif //BUILDER_H
