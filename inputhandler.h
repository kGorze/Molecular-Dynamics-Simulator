#pragma once

#ifndef INPUTHANDLER_H
#define INPUTHANDLER_H

#include <map>
#include <sstream>
#include <vector>
#include <variant>
#include <string>

#include "filehandler.h"


class InputHandler: protected FileHandler {
private:
	std::ifstream mInputDataFile;
	std::ifstream mCoordinatesDataFile;
	std::map<std::string, std::variant<std::string, std::vector<std::string>>> mInputData;
	std::vector<std::vector<std::vector<double>>> mCoordinatesData;

public:
	InputHandler(); //do poprawy
	
	void setVariables(const char* inputFileName, const char* coordinatesFileName);
	void closeFiles();

	void readInputDataFile();
	void loadInputFile(const char* inputFileName);

	void loadCoordinatesFile(const char* coordinatesFileName);
	void readCoordinatesFile();

	void printInputData() const;
	std::map<std::string, std::variant<std::string, std::vector<std::string>>> getmInputData() const;
	std::vector<std::vector<std::vector<double>>> getmCoordinatesData() const;

};

#endif