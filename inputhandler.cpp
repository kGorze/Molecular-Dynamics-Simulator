#include "inputhandler.h"

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


InputHandler::InputHandler() {
}

InputHandler::InputHandler(const char* inputFileName) {

	//could have done it with .good() but this is more readable
	if (!checkFileState(inputFileName)) {
		std::cout << "Error: File " << inputFileName << " not found." << "\n";
		exit(EXIT_FAILURE);
	}

	this->mInputDataFile.open(inputFileName);

	if (!mInputDataFile.is_open()) {
		std::cout << "Error: File " << inputFileName << " could not be opened." << "\n";
		exit(EXIT_FAILURE);
	}

	readInputFile();

	this->mInputDataFile.close();

}
void InputHandler::readInputFile() {
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

	readInputFile();

	this->mInputDataFile.close();
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
