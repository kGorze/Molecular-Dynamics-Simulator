#pragma once
#include <iostream>
#include "filehandler.h"
class OutputHandler: protected FileHandler{
private:
	FILE* mpGeneralDataOutput;
	FILE* mpCoordinatesDataOutput;
	FILE* mpHistogramDataOutput;
public:

};