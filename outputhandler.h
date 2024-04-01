#pragma once
#ifndef OUTPUTHANDLER_H
#define OUTPUTHANDLER_H

#include "filehandler.h"
#include "types.h"

#include <vector>

class OutputHandler: protected FileHandler{
private:
	std::ofstream mGeneralDataOutput;
	std::string mGeneralDataOutputFileName;

	std::ofstream mCoordinatesDataOutput;
	std::string mCoordinatesDataOutputFileName;

	std::ofstream mHistogramDataOutput;
	std::string mHistogramDataOutputFileName;
public:
	
	OutputHandler(); //do poprawy
	~OutputHandler();

	void setVariables(const char* outputFileName, const char* coordinatesFileName, const char* hisogramFileName);
	void closeFiles();

	void writeCoordinatesToFile(const std::vector<std::vector<std::vector<double>>>& data);
	void writeSummaryToFile(int mstepCount, real mtimeNow, VecR mvSum, real mnMol, Prop mtotEnergy, Prop mkinEnergy, Prop mpressure);
	void writeHistogramToFile(const std::vector<real>& histVel);
};

#endif