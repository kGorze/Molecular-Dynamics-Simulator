//
// Created by konrad_guest on 16/05/2024.
//
#include "Headers/director.h"

// Handler for parsing INI file
static int handler(void* user, const char* section, const char* name, const char* value) {
    std::unordered_map<std::string, std::string>* config =
        reinterpret_cast<std::unordered_map<std::string, std::string>*>(user);
    std::string key = std::string(section) + "." + std::string(name);
    (*config)[key] = std::string(value);
    return 1;
}

void Director::set_builder(std::shared_ptr<Builder> b) {
    builder = b;
}

std::shared_ptr<Builder> Director::getBuilder() const {
    return builder;
}

void Director::initializeSimulationWithConfig(const std::string& configFilePath) {
    std::unordered_map<std::string, std::string> config;
    if (ini_parse(configFilePath.c_str(), handler, &config) < 0) {
        std::cout << "Can't load '" << configFilePath << "'\n";
        return;
    }
    try {
        double deltaT = std::stod(config["Parameters.deltaT"]);
        double density = std::stod(config["Parameters.density"]);
        int initUcellX = std::stoi(config["Parameters.initUcell"].substr(0, config["Parameters.initUcell"].find(' ')));
        int initUcellY = std::stoi(config["Parameters.initUcell"].substr(config["Parameters.initUcell"].find(' ') + 1));
        int stepAvg = std::stoi(config["Parameters.stepAvg"]);
        int stepEquil = std::stoi(config["Parameters.stepEquil"]);
        int stepLimit = std::stoi(config["Parameters.stepLimit"]);
        double temperature = std::stod(config["Parameters.temperature"]);
        int limitVel = std::stoi(config["Parameters.limitVel"]);
        double rangeVel = std::stod(config["Parameters.rangeVel"]);
        int sizeHistVel = std::stoi(config["Parameters.sizeHistVel"]);
        int stepVel = std::stoi(config["Parameters.stepVel"]);
        int randSeed = std::stoi(config["Parameters.randSeed"]);
        int numberOfDimensions = std::stoi(config["Parameters.numberOfDimensions"]);

        Eigen::VectorXd initUcell(2);
        initUcell << initUcellX, initUcellY;


        this->builder->setDeltaT(deltaT);
        this->builder->setDensity(density);
        this->builder->setInitUcell(initUcell);
        this->builder->setStepAvg(stepAvg);
        this->builder->setStepEquil(stepEquil);
        this->builder->setStepLimit(stepLimit);
        this->builder->setTemperature(temperature);
        this->builder->setLimitVel(limitVel);
        this->builder->setRangeVel(rangeVel);
        this->builder->setSizeHistVel(sizeHistVel);
        this->builder->setStepVel(stepVel);
        this->builder->setRandSeed(randSeed);
        this->builder->setNumberOfDimensions(numberOfDimensions);
    } catch (const std::out_of_range& e) {
        std::cerr << "Missing parameter in config: " << e.what() << std::endl;
        // Optionally handle missing parameters
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument in config: " << e.what() << std::endl;
        // Optionally handle invalid arguments
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        // Optionally handle unexpected errors
    }
        try {
            this->builder->setConfig(config);
        } catch (const std::exception& e) {
            std::cerr << "Error setting configuration: " << e.what() << std::endl;
        }
}

/* void Director::setRestOfParameters() {

}*/

void Director::buildMinimalSimulation() {
    this->builder->setRestOfParameters();
    this->builder->initializeVectors();
    this->builder->setupStaticSimulation();



    //this->builder->setAtomSize(this->builder->getinitUcellX(), this->builder->getinitUcellY());

    //this->builder->initializeAtoms();
}

void Director::buildFullSimulation() {
    //this->builder->setAtomSize(5,5);
    //this->builder->initializeAtoms();
}

