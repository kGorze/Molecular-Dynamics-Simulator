#pragma once
#include <iostream>
#include "filehandler.h"
class OutputHandler: protected FileHandler{
private:
	std::ifstream mpGeneralDataOutput;
	std::ifstream mpCoordinatesDataOutput;
	std::ifstream mpHistogramDataOutput;
public:
	OutputHandler();
	OutputHandler(const char* outputFileName);

};