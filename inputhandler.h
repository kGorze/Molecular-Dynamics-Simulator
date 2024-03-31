#pragma once
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

#include "filehandler.h"


class InputHandler: protected FileHandler {
private:
	std::ifstream mInputDataFile;
	std::map<std::string, std::variant<std::string, std::vector<std::string>>> mInputData;
public:
	InputHandler(const char* inputFileName);
	void readInputFile();
	void printInputData() const;

};