//WRITTEN BY: github.com/kgorze

#ifndef CONCRETECOMMANDS_H
#define CONCRETECOMMANDS_H


#include "Headers/simulation/simulation2dap.h"
#include "Headers/simulation/builder.h"
#include "Headers/simulation/director.h"
#include "Headers/core/command.h"

class StartSimulationCommand : public command {
private:
    std::shared_ptr<Simulation> simulation;
public:
    StartSimulationCommand(std::shared_ptr<Simulation> sim) : simulation(sim) {}
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
    std::shared_ptr<Simulation> simulation;
public:
    explicit RetrieveDataCommand(std::shared_ptr<Simulation> sim) : simulation(std::move(sim)) {}
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
    std::shared_ptr<Simulation> simulation;
public:
    explicit SaveToCSVCommand(std::shared_ptr<Simulation> sim) : simulation(std::move(sim)) {}
    /* void execute() override {
        if (simulation) {
            /std::ofstream file("simulation_data.csv");
            file << simulation->getAtomsAsCSV();
            file.close();
        }
    }*/
};






#endif //CONCRETECOMMANDS_H
