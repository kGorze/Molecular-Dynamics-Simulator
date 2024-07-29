//
// Created by konrad_guest on 18/05/2024.
//

#ifndef MENU_H
#define MENU_H

#include "Headers/simulation/simulation.h"
#include "Headers/core/command.h"

#include <unordered_map>
#include <memory>
#include <string>
#include <iostream>

class menu {
private:
    std::unordered_map<std::string, std::shared_ptr<command>> commands;
    std::unordered_map<std::string, std::shared_ptr<Simulation>> simulations;

public:
    void setCommand(const std::string& name, std::shared_ptr<command> command);
    void invoke(const std::string& name);
    void setSimulation(const std::string& name, std::shared_ptr<Simulation> simulation);
    std::shared_ptr<command> getCommand(const std::string& name);
    std::shared_ptr<Simulation> getSimulation(const std::string& name);
};



#endif //MENU_H
