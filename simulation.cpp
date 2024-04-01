#pragma once


#include "simulation.h"
#include "rand2D.h"
#include "progressbar.h"

#include <string>
#include <map>
#include <variant>


#include "statistics.h"
#include "leapfrog.h"
#include "physical.h"
#include "vector_operations.h"


#include <chrono>
#include <iostream>
#include <thread>

Simulation::Simulation() {
}


void Simulation::setParameters() {
    std::map<std::string, std::variant<std::string, std::vector<std::string>>> inputMap = minputHandler.getmInputData();

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

void Simulation::initVectors() {
    DO_MOL{
        this->mMolsVector.push_back(Mol());
    }
    for(int i = 0;i< this->msizeHistVel;i++){
		this->mhistogramVelocities.push_back(0);
	}
    for(int i = 0;i< this->mstepLimit;i++){
		this->mdataCoordinates.push_back(std::vector<std::vector<double> >());
	}
}

void Simulation::initCoordinates() {
    VecR c, gap;
    gap.x = 0;
    gap.y = 0;

    c.x = 0;
    c.y = 0;

    int nx, ny;
    VDiv(gap, mregion, minitUcell);
    this->mnDOMOL = 0;
    for (nx = 0; nx < minitUcell.x; nx++) {
        for (ny = 0; ny < minitUcell.y; ny++) {
			VSet(c, nx + 0.5, ny + 0.5);
			VMul(c, c, gap);
			VVSAdd(c, -0.5, mregion); 
			this->mMolsVector.at(this->mnDOMOL).coordinates = c;
            ++this->mnDOMOL;
        }
    }
}

void Simulation::initVelocities() {
    VZero(this->mvSum);
    DO_MOL{
        velocity_rand(this->mMolsVector.at(this->mnDOMOL));
    VScale(this->mMolsVector.at(this->mnDOMOL).velocity, mvelMag);
    VVAdd(this->mvSum, this->mMolsVector.at(this->mnDOMOL).velocity);
    }
        DO_MOL{
            VVSAdd(this->mMolsVector.at(this->mnDOMOL).velocity, (-1.0 / this->mnMol), this->mvSum);
    }
};

void Simulation::initAccelerations() {
    DO_MOL{
        VZero(this->mMolsVector.at(this->mnDOMOL).accelaration);
    }
};

void Simulation::accumulateProperties(unsigned int icode) {
    if (icode == 0) {
        PropZero(this->mtotEnergy);
        PropZero(this->mkinEnergy);
        PropZero(this->mpressure);
    }
    else if (icode == 1) {
        PropAccum(this->mtotEnergy);
        PropAccum(this->mkinEnergy);
        PropAccum(this->mpressure);
    }
    else if (icode == 2) {
        PropAvg(this->mtotEnergy, this->mstepAvg);
        PropAvg(this->mkinEnergy, this->mstepAvg);
        PropAvg(this->mpressure, this->mstepAvg);
    }
};

void Simulation::setupJob() {
    this->mstepCount = 0;
    this->mcountVel = 0;
    this->initCoordinates();
    this->initVelocities();
    this->initAccelerations();
    this->accumulateProperties(0);
};

void Simulation::leapFrogStep(unsigned int part) {
    if (part == 1) {
        DO_MOL{
            leapfrog_velocity(this->mMolsVector.at(this->mnDOMOL).velocity, 0.5 * this->mdeltaT, this->mMolsVector.at(this->mnDOMOL).accelaration);
            leapfrog_coordinates(this->mMolsVector.at(this->mnDOMOL).coordinates, this->mdeltaT, this->mMolsVector.at(this->mnDOMOL).velocity);
        }
    }
    else {
        DO_MOL{
            leapfrog_velocity(this->mMolsVector.at(this->mnDOMOL).velocity, 0.5 * this->mdeltaT, this->mMolsVector.at(this->mnDOMOL).accelaration);
        }
    }
}

void Simulation::applyBoundaryConditions() {
        DO_MOL{
            VWrapAll(this->mMolsVector.at(this->mnDOMOL).coordinates, this->mregion);
        }
}

void Simulation::computeForces() {
        VecR dr;
        dr.x = 0;
        dr.y = 0;

        //fcVal - Lennard-Jones force
        //rr - distance between two molecules
        //rrCut - cut-off distance
        //rri - 1/rr
        //rri3 - 1/rr^3

        real fcVal, rr, rrCut, rri, rri3; 
        int j1, j2, n; 

        rrCut = Sqr(this->mrCut);

        DO_MOL{
            this->mMolsVector.at(this->mnDOMOL).accelaration.x = 0.;
            this->mMolsVector.at(this->mnDOMOL).accelaration.y = 0.;
        }

        this->muSum = 0.;
        this->mvirSum = 0;

        for (j1 = 0; j1 < this->mnMol - 1; j1++) {
            for (j2 = j1 + 1; j2 < this->mnMol; j2++) {

                dr.x = this->mMolsVector.at(j1).coordinates.x - this->mMolsVector.at(j2).coordinates.x;
                dr.y = this->mMolsVector.at(j1).coordinates.y - this->mMolsVector.at(j2).coordinates.y;

                if (dr.x >= 0.5 * this->mregion.x) dr.x -= this->mregion.x;
                else if (dr.x < -0.5 * this->mregion.x) dr.x += this->mregion.x;
                if (dr.y >= 0.5 * this->mregion.y) dr.y -= this->mregion.y;
                else if (dr.y < -0.5 * this->mregion.y) dr.y += this->mregion.y;

                rr = dr.x * dr.x + dr.y * dr.y;

                if (rr < rrCut) {
                    rri = 1. / rr;
                    rri3 = rri * rri * rri;

                    fcVal = 48. * rri3 * (rri3 - 0.5) * rri;


                    mMolsVector.at(j1).accelaration.x += fcVal * dr.x;
                    mMolsVector.at(j1).accelaration.y += fcVal * dr.y;
                    mMolsVector.at(j2).accelaration.x -= fcVal * dr.x;
                    mMolsVector.at(j2).accelaration.y -= fcVal * dr.y;
                    
                    this->muSum += 4. * rri3 * (rri3 - 1.) + 1.;
                    this->mvirSum += fcVal * rr;
                }
            }
        }
};

void Simulation::evaluateProperties() {
        real vv; 

        VZero(this->mvSum);

        this->mvvSum = 0;

        DO_MOL{
            VVAdd(this->mvSum, this->mMolsVector.at(this->mnDOMOL).velocity);
            vv = VLenSq(this->mMolsVector.at(this->mnDOMOL).velocity);
            this->mvvSum += vv;
        }
        this->mkinEnergy.val = (0.5 * this->mvvSum) / this->mnMol;
        this->mtotEnergy.val = this->mkinEnergy.val + (this->muSum / this->mnMol);
        this->mpressure.val = this->mdensity * (this->mvvSum - this->mvirSum) / (this->mnMol * NDIM);



};

void Simulation::evaluateVelocityDistribution() {
        this->mdeltaV = 0;
        this->mhistSum = 0;
        int j;

        if (this->mcountVel == 0) {
            for (j = 0; j < this->msizeHistVel; j++) {
                this->mhistogramVelocities.at(j) = 0;
            }
        }
        this->mdeltaV = this->mrangeVel / this->msizeHistVel;
        DO_MOL{
            j = (VLen(this->mMolsVector.at(this->mnDOMOL).velocity)) / this->mdeltaV;
            this->mhistogramVelocities.at(Min(j, this->msizeHistVel - 1)) += 1;
        }
        this->mcountVel += 1;
        if (this->mcountVel == this->mlimitVel) {
            this->mhistSum = 0;
            for (j = 0; j < this->msizeHistVel; j++) {
                this->mhistSum += this->mhistogramVelocities.at(j);
            }
            for (j = 0; j < this->msizeHistVel; j++) {
                this->mhistogramVelocities.at(j) /= this->mhistSum;
            }
            this->mcountVel = 0;
        }
        this->mhFunction = 0;
        for (j = 0; j < this->msizeHistVel; j++) {
            if (this->mhistogramVelocities.at(j) > 0) {
                this->mhFunction += this->mhistogramVelocities.at(j) * log(this->mhistogramVelocities.at(j) / ((j + 0.5) * this->mdeltaV));
            }
        }
};

void Simulation::printSummaryToScreen() {

        real vSumValue = VCSum(this->mvSum) / this->mnMol;
        double totEnergyEst = PropEst(this->mtotEnergy);
        double totEnergySig = PropEstSig(this->mtotEnergy);
        double kinEnergyEst = PropEst(this->mkinEnergy);
        double kinEnergySig = PropEstSig(this->mkinEnergy);
        double pressureEst = PropEst(this->mpressure);
        double pressureSig = PropEstSig(this->mpressure);

        std::cout << "\n-----------------------------------------------------------------------\n";
        std::cout << "Simulation result\n";
        std::cout << "-----------------------------------------------------------------------\n\n";

        std::cout << " Step   Time    Sum(v)  Etot            EtotSig            Ekin            EkinSig            Pressure        PressureSig\n";

        std::cout << this->mstepCount << " " << this->mtimeNow << " " << (VCSum(this->mvSum) / this->mnMol) << " "
            << totEnergyEst << " " << totEnergySig
            << " " << kinEnergyEst << " " << kinEnergySig
            << " " << pressureEst << " " << pressureSig << "\n";
 };

void Simulation::saveSummaryToFile() {
	this->moutputHandler.writeSummaryToFile(this->mstepCount, this->mtimeNow, this->mvSum, this->mnMol, this->mtotEnergy, this->mkinEnergy, this->mpressure);
}

double Simulation::singleSimulationStep() {
        this->mstepCount++;
        this->mtimeNow = this->mstepCount * this->mdeltaT;

        clock_t start = clock();


        this->leapFrogStep(1);
        this->applyBoundaryConditions();
        this->computeForces();
        this->leapFrogStep(2);

        this->evaluateProperties();
        this->accumulateProperties(1);
        DO_MOL{
            this->mdataCoordinates.at(this->mstepCount - 1).push_back({ this->mMolsVector.at(this->mnDOMOL).coordinates.x, this->mMolsVector.at(this->mnDOMOL).coordinates.y});
        }
        if (this->mstepCount % this->mstepAvg == 0) {
            this->accumulateProperties(2);
            if (this->mstepCount >= this->mstepEquil && (this->mstepCount - this->mstepEquil) % this->mstepVel == 0) {
                this->evaluateVelocityDistribution();
            }
            //PrintSummary(fp);
            if(this->mstepCount != this->mstepLimit){ this->accumulateProperties(0); }

        }
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        return time_taken;

 };


std::vector<std::vector<std::vector<double>>> Simulation::readInputCoordiantesFile(const char* coordinatesFileName) {
    this->minputHandler.loadCoordinatesFile(coordinatesFileName);
    return minputHandler.getmCoordinatesData();
};

void Simulation::saveCoordinatesToFile() {
    moutputHandler.writeCoordinatesToFile(this->mdataCoordinates);
}

void Simulation::saveHistogramToFile() {
	moutputHandler.writeHistogramToFile(this->mhistogramVelocities);
}

void Simulation::runSimulation(unsigned int option) {
        std::cout << "-----------------------------------------------------------------------\n";
        std::cout << "Initializing simulation\n";
        std::cout << "-----------------------------------------------------------------------\n\n";
        minputHandler.loadInputFile("data.in");

        std::cout << "-----------------------------------------------------------------------\n";
        std::cout << "Parameters\n";
        std::cout << "-----------------------------------------------------------------------\n";
        minputHandler.printInputData();
        std::cout << "-----------------------------------------------------------------------\n\n";

        //std::this_thread::sleep_for(std::chrono::seconds(4));
        
        this->setParameters();

        this->initVectors();
        this->initCoordinates();
        this->setupJob();

        system("cls");

        uint64_t total_steps = mstepLimit;
        this->mprogressBar.totalTicks(total_steps);

        int moreCycles = 1;

        while (moreCycles) {
            this->singleSimulationStep();
            if (this->mstepCount >= this->mstepLimit) {
                moreCycles = 0;
            }
            this->mprogressBar.tick();
        }
        system("cls");

        if (option == 1) {
            this->printSummaryToScreen();
            this->saveSummaryToFile();
            this->saveCoordinatesToFile();
            this->saveHistogramToFile();
        }
        minputHandler.closeFiles();
        moutputHandler.closeFiles();
};

void Simulation::setFiles(const char* inputFileName, const char* coordinatesFileName, const char* outputFileName, const char* histogramFileName) {
	this->minputHandler.setVariables(inputFileName, coordinatesFileName);
	this->moutputHandler.setVariables(outputFileName, coordinatesFileName, histogramFileName);
}

std::vector<std::vector<std::vector<double>>> Simulation::getCoordinatesData() {
    return this->mdataCoordinates;
}

void Simulation::setCoordinatesData(std::vector<std::vector<std::vector<double>>> &dataCoords) {
    this->mdataCoordinates = dataCoords;
}