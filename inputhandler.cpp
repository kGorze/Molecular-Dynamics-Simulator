#pragma once

#include "inputhandler.h"
#include "filehandler.h"


#include <vector>
#include <sstream>
#include <variant>
#include <string>

InputHandler::InputHandler() {
}

void InputHandler::setVariables(const char* inputFileName, const char* coordinatesFileName) {
	if (!checkFileState(inputFileName)) {
		std::cout << "Error: File " << inputFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	this->mInputDataFile.open(inputFileName);

	if (!mInputDataFile.is_open()) {
		std::cout << "Error: File " << inputFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	readInputDataFile();

	if (!checkFileState(coordinatesFileName)) {
		std::cout << "Error: File " << coordinatesFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}
	this->mCoordinatesDataFile.open(coordinatesFileName);

	if (!mCoordinatesDataFile.is_open()) {
		std::cout << "Error: File " << coordinatesFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}
	//INICJALIZACJA WEKTORA KTÓRY JEST PUSTY - DO POPRAWY

	//readCoordinatesFile();
	//std::vector<std::vector<std::vector<double>>> temp = this->getmCoordinatesData();
}

void InputHandler::readInputDataFile() {
	std::string line;
	std::string key, value;
	std::vector<std::string> values;

	while (std::getline(mInputDataFile, line)) {
		std::istringstream iss(line);
		iss >> key;


		while (iss >> value) {
			values.push_back(value);
		}

		if (values.size() == 1) this->mInputData[key] = values[0];
		else this->mInputData[key] = values;


		key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
		for (auto& value : values)
			value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());
		values.clear();
	}
}

void InputHandler::loadInputFile(const char* inputFileName) {

	if (!checkFileState(inputFileName)) {
		std::cout << "Error: File " << inputFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	this->mInputDataFile.open(inputFileName);

	if (!mInputDataFile.is_open()) {
		std::cout << "Error: File " << inputFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	readInputDataFile();

	this->mInputDataFile.close();
}


void InputHandler::readCoordinatesFile() {
	int iteration;
	double x, y;

	while (mCoordinatesDataFile >> iteration >> x >> y) {
		if (mCoordinatesDataFile.fail()) {
			std::cerr << "Error: Failed to read data from file." << std::endl;
			mCoordinatesDataFile.clear();
			mCoordinatesDataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		mCoordinatesData.at(iteration).push_back({ x, y });
	}
}

void InputHandler::loadCoordinatesFile(const char* coordinatesFileName) {
	if (!checkFileState(coordinatesFileName)) {
		std::cout << "Error: File " << coordinatesFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}
	this->mCoordinatesDataFile.open(coordinatesFileName);
	if (!mCoordinatesDataFile.is_open()) {
		std::cout << "Error: File " << coordinatesFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}
	readCoordinatesFile();
	this->mCoordinatesDataFile.close();
}




void InputHandler::printInputData() const {
	for (auto& [key, entry] : this->mInputData) {
		if (std::holds_alternative<std::string>(entry)) {
			std::cout << "Key: " << key << "\n";
			std::cout << "Value: " << std::get<std::string>(entry) << "\n";
		}
		else {
			std::cout << "Key: " << key << "\n";
			std::cout << "Values:" << "\n";
			for (const auto& val : std::get<std::vector<std::string>>(entry)) {
				std::cout << val << "\t";
			}
			std::cout << "\n";
		}
	}
}
std::map<std::string, std::variant<std::string, std::vector<std::string>>> InputHandler::getmInputData() const {
	return this->mInputData;
}
std::vector<std::vector<std::vector<double>>> InputHandler::getmCoordinatesData() const {
	return this->mCoordinatesData;
}

void InputHandler::closeFiles() {
	if (this->mInputDataFile.is_open()) {
		this->mInputDataFile.close();
	}
	if (this->mCoordinatesDataFile.is_open()) {
		this->mCoordinatesDataFile.close();
	}
}