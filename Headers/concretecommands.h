//
// Created by konrad_guest on 18/05/2024.
//

#ifndef CONCRETECOMMANDS_H
#define CONCRETECOMMANDS_H



#include "Command.h"
#include "builder.h"
#include "director.h"

class StartSimulationCommand : public command {
private:
    std::shared_ptr<Simulation2D> simulation;
public:
    StartSimulationCommand(std::shared_ptr<Simulation2D> sim) : simulation(sim) {}
    /*void execute() override {
        simulation->run();
    }*/
};

class InitializeSimulationCommand : public command {
private:
    Director& director;
public:
    InitializeSimulationCommand(Director& dir) : director(dir) {}
    void execute() override {
        director.buildMinimalSimulation();
    }
};

class RetrieveDataCommand : public command {
private:
    std::shared_ptr<Simulation2D> simulation;
public:
    explicit RetrieveDataCommand(std::shared_ptr<Simulation2D> sim) : simulation(std::move(sim)) {}
    void execute() override {
        if (simulation) {
            //simulation->printAtoms();
        }
    }
    /* std::string getResults() const {
        return simulation ? simulation->getAtomsAsString() : "";
    }*/
};

class SaveToCSVCommand : public command {
private:
    std::shared_ptr<Simulation2D> simulation;
public:
    explicit SaveToCSVCommand(std::shared_ptr<Simulation2D> sim) : simulation(std::move(sim)) {}
    /* void execute() override {
        if (simulation) {
            /std::ofstream file("simulation_data.csv");
            file << simulation->getAtomsAsCSV();
            file.close();
        }
    }*/
};






#endif //CONCRETECOMMANDS_H
