#pragma once
#include <iostream>
#include "types.h"
#include "OutputHandler.h"
#include "InputHandler.h"
#include "simulation.h"
#include "vector_operations.h"
#include "rand2D.h"

Simulation::Simulation() {
}

Simulation::Simulation(const char* outputFileName, const char* inputFileName) {
    //moutputHandler = OutputHandler(outputFileName);
    this->mkeyValueInputHandler.loadInputFile(inputFileName);
}

void Simulation::setParameters() {
    std::map<std::string, std::variant<std::string, std::vector<std::string>>> inputMap = mkeyValueInputHandler.getmInputData();

    //const char* variableNames[] = {
    //    "mregion", "mvSum", "minitUcell","mrCut", "mkinEnergy", "mtotEnergy", "mpressure", "mdeltaT", "mdensity", "mtemperature", "mtimeNow",
    //    "muCut", "muSum", "mvelMag", "mvirSum", "mvvSum", "mrangeVel", "mhFunction", "mmoreCycles", "mstepAvg", "mstepCount",
    //    "mstepEquil", "mstepLimit", "mcountVel", "mlimitVel", "mstepVel", "mnMol", "msizeHistVel"
    //};

    if (inputMap.find("deltaT") != inputMap.end()) {
        mdeltaT = std::stod(std::get<std::string>(inputMap["deltaT"]));
    }
    if (inputMap.find("density") != inputMap.end()) {
		mdensity = std::stod(std::get<std::string>(inputMap["density"]));
	}
    if (inputMap.find("initUcell") != inputMap.end()) {
		std::vector<std::string> initUcell = std::get<std::vector<std::string>>(inputMap["initUcell"]);
		minitUcell.x = std::stod(initUcell[0]);
        minitUcell.y = std::stod(initUcell[1]);
	}
    if (inputMap.find("stepAvg") != inputMap.end()) {
		mstepAvg = std::stoi(std::get<std::string>(inputMap["stepAvg"]));
	}
    if (inputMap.find("stepEquil") != inputMap.end()) {
		mstepEquil = std::stoi(std::get<std::string>(inputMap["stepEquil"]));
	}
    if (inputMap.find("stepLimit") != inputMap.end()) {
		mstepLimit = std::stoi(std::get<std::string>(inputMap["stepLimit"]));
	}
    if (inputMap.find("temperature") != inputMap.end()) {
		mtemperature = std::stod(std::get<std::string>(inputMap["temperature"]));
	}
    if (inputMap.find("limitVel") != inputMap.end()) {
		mlimitVel = std::stoi(std::get<std::string>(inputMap["limitVel"]));
	}
    if (inputMap.find("rangeVel") != inputMap.end()) {
		mrangeVel = std::stod(std::get<std::string>(inputMap["rangeVel"]));
	}
    if (inputMap.find("sizeHistVel") != inputMap.end()) {
		msizeHistVel = std::stoi(std::get<std::string>(inputMap["sizeHistVel"]));
	}
    if (inputMap.find("stepVel") != inputMap.end()) {
		mstepVel = std::stoi(std::get<std::string>(inputMap["stepVel"]));
	}
    if (inputMap.find("randSeed") != inputMap.end()) {
		mrandSeed = std::stoi(std::get<std::string>(inputMap["randSeed"]));
	}

    this->mrCut = pow(2.0, 1.0 / 6.0);
    VSCopy(this->mregion, 1 / sqrt(this->mdensity), this->minitUcell);
    this->mnMol = VProd(this->minitUcell);
    this->mvelMag = sqrt(NDIM * (1 - (1 / this->mnMol) * this->mtemperature));
}

void Simulation::initCoordinates() {
    VecR c, gap;
    gap.x = 0;
    gap.y = 0;

    c.x = 0;
    c.y = 0;

    int n, nx, ny;
    VDiv(gap, mregion, minitUcell);
    n = 0;
    for (nx = 0; nx < minitUcell.x; nx++) {
        for (ny = 0; ny < minitUcell.y; ny++) {
			VSet(c, nx + 0.5, ny + 0.5);
			VMul(c, c, gap);
			VVSAdd(c, -0.5, mregion); 
			this->mol.at(n).coordinates = c;
			++n;
		}
	}
}

void Simulation::initVelocities() {
        int n;
        VZero(this->mvSum);
        DO_MOL{
            velocity_rand(this->mol.at(n));
        VScale(this->mol.at(n).velocity, mvelMag);
        VVAdd(this->mvSum, this->mol.at(n).velocity);
        }
            DO_MOL{
                VVSAdd(this->mol.at(n).velocity, (-1.0 / this->mnMol), this->mvSum);
        }
    };

void Simulation::setupJob() {
    this->mstepCount = 0;
    this->mcountVel = 0;
    this->initCoordinates(); // initial coordinates
    this->initVelocities(); // initial velocities

    //ZROBIĆ  N-MOL

    InitAccels(); // initial accelerations
    AccumProps(0);
};