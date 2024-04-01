#pragma once
#include <iostream>
#include "types.h"
#include "OutputHandler.h"
#include "InputHandler.h"
#include "vector_operations.h"
class Simulation {
private:

	OutputHandler moutputHandler;
	InputHandler mkeyValueInputHandler;


	VecR mregion, mvSum, minitUcell;
	Prop mkinEnergy, mtotEnergy, mpressure;
	real mdeltaT, mdensity,mrCut, mtemperature, mtimeNow, muCut, muSum, mvelMag, mvirSum, mvvSum, mrangeVel, mhFunction, mrandSeed;
	int mmoreCycles, mstepAvg, mstepCount, mstepEquil, mstepLimit, mcountVel, mlimitVel, mstepVel, mnMol, msizeHistVel;

	//std::vector<Mol> mol(nMol);
	std::vector<Mol> mol;
	//std::vector<real> histVel(sizeHistVel);
	std::vector<real> histVel;

public:
	Simulation();
	Simulation(const char* outputFileName, const char* inputFileName);
	void setParameters();
	void initCoordinates();
	void initVelocities();
	void setupJob();


};