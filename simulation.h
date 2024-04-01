#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include "outputhandler.h"
#include "inputhandler.h"
#include "vector_operations.h"

#include <vector>

class Simulation {
private:

	OutputHandler moutputHandler;
	InputHandler minputHandler;


	VecR mregion, mvSum, minitUcell;
	Prop mkinEnergy, mtotEnergy, mpressure;
	real mdeltaT, mdensity,mrCut, mtemperature, mtimeNow, muCut, muSum, mvelMag, mvirSum, mvvSum, mrangeVel, mhFunction, mrandSeed,mdeltaV,mhistSum;
	int mmoreCycles, mstepAvg, mstepCount, mstepEquil, mstepLimit, mcountVel, mlimitVel, mstepVel, mnMol, msizeHistVel,mnDOMOL;

	//std::vector<Mol> mol(nMol);
	std::vector<Mol> mMolsVector;
	//std::vector<real> histVel(sizeHistVel);
	std::vector<real> mhistogramVelocities;
	//std::vector<std::vector<std::vector<double>>> dataCoords(500);
	std::vector<std::vector<std::vector<double>>> mdataCoordinates;

public:
	Simulation();

	void setFiles(const char* inputFileName, const char* coordinatesFileName, const char* outputFileName, const char* histogramFileName);

	void printSummaryToScreen();
	void saveSummaryToFile();

	void saveCoordinatesToFile();

	void saveHistogramToFile();


	std::vector<std::vector<std::vector<double>>> readInputCoordiantesFile(const char* coordinatesFileName);


	void setParameters();

	void initVectors();
	void initCoordinates();
	void initVelocities();
	void initAccelerations();
	void accumulateProperties(unsigned int icode);
	void setupJob();
	void leapFrogStep(unsigned int part);
	void applyBoundaryConditions();
	void computeForces();
	void evaluateProperties();
	void evaluateVelocityDistribution();
	double singleSimulationStep();

	std::vector<std::vector<std::vector<double>>> getCoordinatesData();

	void runSimulation(unsigned int option);


};

#endif	