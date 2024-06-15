//
// Created by konrad_guest on 18/05/2024.
//

#include "Headers/menu.h"

void menu::setCommand(const std::string& name, std::shared_ptr<command> command) {
    commands[name] = command;
}

void menu::invoke(const std::string& name) {
    if (commands.find(name) != commands.end()) {
        commands[name]->execute();
    }
}

void menu::setSimulation(const std::string& name, std::shared_ptr<Simulation> simulation) {
    simulations[name] = simulation;
}

std::shared_ptr<command> menu::getCommand(const std::string& name) {
    auto it = commands.find(name);
    if (it != commands.end()) {
        return it->second;
    }
    return nullptr;

}

std::shared_ptr<Simulation> menu::getSimulation(const std::string& name) {
    if (simulations.find(name) != simulations.end()) {
        return simulations[name];
    }
    return nullptr;
}