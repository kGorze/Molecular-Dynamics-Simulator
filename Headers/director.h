//
// Created by konrad_guest on 16/05/2024.
//

#ifndef DIRECTOR_H
#define DIRECTOR_H
#include "builder.h"
#include "Third_party/inih/ini.h"
#include <memory>
#include <unordered_map>

class Director {
private:
    std::shared_ptr<Builder> builder;
public:
    void set_builder(std::shared_ptr<Builder> b);

    void initializeSimulationWithConfig(const std::string& configFilePath);
    void buildMinimalSimulation();
    void buildFullSimulation();

    std::shared_ptr<Builder> getBuilder() const;
};




#endif //DIRECTOR_H
