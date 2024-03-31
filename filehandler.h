#pragma once
#include <iostream>
#include <fstream>

class FileHandler {
protected:
    bool checkFileState(const char* fileName) {
        std::ifstream file(fileName);
        return file.good();
    }
};