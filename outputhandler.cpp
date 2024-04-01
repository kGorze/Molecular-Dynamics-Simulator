#pragma once

#include "outputhandler.h"
#include "vector_operations.h"
#include "statistics.h"

#include <vector>
#include <iomanip>

OutputHandler::OutputHandler() {
}

void OutputHandler::setVariables(const char* outputFileName, const char* coordinatesFileName, const char* hisogramFileName) {
	if (!checkFileState(outputFileName)) {
		std::cout << "Error: File " << outputFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	mGeneralDataOutput.open(outputFileName);

	if (!mGeneralDataOutput.is_open()) {
		std::cout << "Error: File " << outputFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	this->mGeneralDataOutputFileName = outputFileName;

	if (!checkFileState(coordinatesFileName)) {
		std::cout << "Error: File " << coordinatesFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	mCoordinatesDataOutput.open(coordinatesFileName);

	if (!mCoordinatesDataOutput.is_open()) {
		std::cout << "Error: File " << coordinatesFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	this->mCoordinatesDataOutputFileName = coordinatesFileName;

	if (!checkFileState(hisogramFileName)) {
		std::cout << "Error: File " << hisogramFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	mHistogramDataOutput.open(hisogramFileName);

	if (!mHistogramDataOutput.is_open()) {
		std::cout << "Error: File " << hisogramFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	this->mHistogramDataOutputFileName = hisogramFileName;

	//do stuff
}

OutputHandler::~OutputHandler(){
	if (this->mGeneralDataOutput.is_open()) {
		this->mGeneralDataOutput.close();
	}
	if (this->mCoordinatesDataOutput.is_open()) {
		this->mCoordinatesDataOutput.close();
	}
	if (this->mHistogramDataOutput.is_open()) {
		this->mHistogramDataOutput.close();
	}

}

void OutputHandler::writeCoordinatesToFile(const std::vector<std::vector<std::vector<double>>>& data) {
	if (!this->mCoordinatesDataOutput.is_open()) {
		std::cerr << "Error: File is not open for writing." << "\n";
		return;
	}

	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data.at(i).size(); j++) {
			this->mCoordinatesDataOutput << std::setw(5) << i <<"\t" << std::setw(8) << std::fixed << std::setprecision(4)
				<< data.at(i).at(j).at(0) << "\t" << std::setw(8) << std::fixed << std::setprecision(4)
				<< data.at(i).at(j).at(1) << "\n";
		}
	}
}

void OutputHandler::writeSummaryToFile(int mstepCount, real mtimeNow, VecR mvSum, real mnMol, Prop mtotEnergy, Prop mkinEnergy, Prop mpressure) {
	if (!this->mGeneralDataOutput.is_open()) {
		std::cerr << "Error: File is not open for writing." << "\n";
		return;
	}

	real vSumValue = VCSum(mvSum) / mnMol;
	real totEnergyEst = PropEst(mtotEnergy);
	real totEnergySig = PropEstSig(mtotEnergy);
	real kinEnergyEst = PropEst(mkinEnergy);
	real kinEnergySig = PropEstSig(mkinEnergy);
	real pressureEst = PropEst(mpressure);
	real pressureSig = PropEstSig(mpressure);

	mGeneralDataOutput << " Step   Time    Sum(v)  Etot            EtotSig            Ekin            EkinSig            Pressure        PressureSig\n";

	mGeneralDataOutput << std::setw(5) << mstepCount << std::setw(8) << std::fixed << std::setprecision(4)
		<< mtimeNow << std::setw(7) << std::fixed << std::setprecision(20)
		<< vSumValue << std::setw(7) << std::fixed << std::setprecision(20)
		<< totEnergyEst << std::setw(7) << std::fixed << std::setprecision(20)
		<< totEnergySig << std::setw(7) << std::fixed << std::setprecision(20)
		<< kinEnergyEst << std::setw(7) << std::fixed << std::setprecision(20)
		<< kinEnergySig << std::setw(7) << std::fixed << std::setprecision(20)
		<< pressureEst << std::setw(7) << std::fixed << std::setprecision(20)
		<< pressureSig << "\n";
	};


void OutputHandler::writeHistogramToFile(const std::vector<real>& histVel) {
	if (!this->mHistogramDataOutput.is_open()) {
		std::cerr << "Error: File is not open for writing." << "\n";
		return;
	}

	for (int i = 0; i < histVel.size(); i++) {
		this->mHistogramDataOutput << std::setw(5) << i <<"\t" << std::setw(8) << std::fixed << std::setprecision(4)
			<< histVel.at(i) << "\n";
	}
}

void OutputHandler::closeFiles() {
	if (this->mGeneralDataOutput.is_open()) {
		this->mGeneralDataOutput.close();
	}
	if (this->mCoordinatesDataOutput.is_open()) {
		this->mCoordinatesDataOutput.close();
	}
	if (this->mHistogramDataOutput.is_open()) {
		this->mHistogramDataOutput.close();
	}
}