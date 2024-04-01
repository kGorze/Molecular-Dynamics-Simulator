#pragma once
#include "filehandler.h"
#include "outputhandler.h"

#ifndef IOSTREAM_H
#define IOSTREAM_H
#include <iostream>
#endif

#ifndef FSTREAM_H
#define FSTREAM_H
#include <fstream>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

#ifndef VARIANT_H
#define VARIANT_H
#include <variant>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef SSTREAM_H
#define SSTREAM_H
#include <sstream>
#endif

OutputHandler::OutputHandler() {
}


OutputHandler::OutputHandler(const char* outputFileName) {
	if (!checkFileState(outputFileName)) {
		std::cout << "Error: File " << outputFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	mpGeneralDataOutput.open(outputFileName);

	if (!mpGeneralDataOutput.is_open()) {
		std::cout << "Error: File " << outputFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	//do stuff

	mpGeneralDataOutput.close();
}
