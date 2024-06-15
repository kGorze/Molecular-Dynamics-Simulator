#include "Headers/director.h"
#include "Headers/builder.h"
#include "Headers/command.h"
#include "Headers/concretecommands.h"
#include "Headers/menu.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include "Third_party/zmq/zmq.hpp"
#include <thread>
#include <atomic>


/*
void listenForCommands(zmq::socket_t& socket, menu& menu, std::atomic<bool>& running) {
    while (running) {
        zmq::message_t reply;
        if (socket.recv(reply, zmq::recv_flags::none)) {
            std::string replyStr(static_cast<char*>(reply.data()), reply.size());
            std::cout << "Received command from server: " << replyStr << std::endl;

            if (replyStr == "initialize") {
                menu.invoke("initialize");
            } else if (replyStr == "start") {
                menu.invoke("start");
            } else if (replyStr == "retrieve") {
                menu.invoke("retrieve");
                auto command = std::dynamic_pointer_cast<RetrieveDataCommand>(menu.getCommand("retrieve"));
                if (command) {
                    std::string data = command->getResults();
                    zmq::message_t resultsMessage(data.size());
                    memcpy(resultsMessage.data(), data.c_str(), data.size());
                    socket.send(resultsMessage, zmq::send_flags::none);
                }
            } else if (replyStr == "save") {
                menu.invoke("save");
                zmq::message_t ackMessage(3);
                memcpy(ackMessage.data(), "ack", 3);
                socket.send(ackMessage, zmq::send_flags::none);
            }

            std::cout << "Command executed and results sent to server" << std::endl;
        }
    }
}
 */

int main() {
    Director director;
    auto builder2D = std::make_shared<Simulation2DBuilder>();
    director.set_builder(builder2D);

    //CODE FOR PATH
    std::string configFilePath = R"(C:\Users\konrad_guest\CLionProjects\MDS\Parameters\config.ini)";
    auto simulation2D = std::dynamic_pointer_cast<Simulation2D>(builder2D->getSimulation());
    if (simulation2D) {
        director.initializeSimulationWithConfig(configFilePath);
        director.buildMinimalSimulation();
        director.getBuilder()->getSimulation()->printConfig();
        director.getBuilder()->getSimulation()->run(1);
        director.getBuilder()->openCoordinatesFile("Resources/results.csv", 1);
        director.getBuilder()->writeCoordinatesToFile(simulation2D->getDataCoordinates());
        director.getBuilder()->closeCoordinatesFile();


        //simulation2D->printAtomsSize();
        //simulation2D->printAtomsSize();
        //simulation2D->printAtoms();
    }

    return 0;
    /*auto startCommand = std::make_shared<StartSimulationCommand>(simulation2D);
    auto initCommand = std::make_shared<InitializeSimulationCommand>(director);
    auto retrieveCommand = std::make_shared<RetrieveDataCommand>(simulation2D);
    auto saveCommand = std::make_shared<SaveToCSVCommand>(simulation2D);*/

    /*menu menu;
    menu.setCommand("start", startCommand);
    menu.setCommand("initialize", initCommand);
    menu.setCommand("retrieve", retrieveCommand);
    menu.setCommand("save", saveCommand);
    menu.setSimulation("Simulation2D", simulation2D);*/

    /*zmq::context_t context(1);
    zmq::socket_t socket(context, zmq::socket_type::dealer);

    std::string identity = "Client" + std::to_string(rand());
    socket.set(zmq::sockopt::identity, identity);

    try {
    socket.connect("tcp://localhost:5555");
    std::cout << "Connected to server with identity: " << identity << std::endl;

    std::atomic<bool> running(true);
    std::thread listenerThread(listenForCommands, std::ref(socket), std::ref(menu), std::ref(running));

    std::cout << "Press Enter to stop..." << std::endl;
    std::cin.get();

    running = false;
    listenerThread.join();

    } catch (const zmq::error_t &e) {
    std::cerr << "ZeroMQ error: " << e.what() << std::endl;
    } catch (const std::exception &e) {
    std::cerr << "Standard exception: " << e.what() << std::endl;
    } catch (...) {
    std::cerr << "Unknown error occurred" << std::endl;
    }*/


}