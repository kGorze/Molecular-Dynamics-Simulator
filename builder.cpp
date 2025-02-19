//WRITTEN BY: github.com/kgorze

#include "Headers/simulation/builder.h"

#include <iostream>


Eigen::Vector2d velocity_rand(std::mt19937& rng, std::uniform_real_distribution<double>& dist) {
    Eigen::Vector2d velocity;
    velocity.x() = std::cos(2 * M_PI * dist(rng));
    velocity.y() = std::cos(2 * M_PI * dist(rng));
    return velocity;
}

Builder::Builder() {
    std::cout << "Builder constructor" << std::endl;
}

Builder::~Builder() {
    std::cout << "Builder destructor" << std::endl;
}


void Simulation2DBuilder::ensureDirectoryExists(const std::string &path) {
    std::filesystem::path contentRoot = std::filesystem::current_path();
    std::filesystem::path dirPath = std::filesystem::path(path).parent_path();
    if (!std::filesystem::exists(dirPath)) {
        std::filesystem::create_directories(dirPath);
    }
}


void Simulation2DBuilder::setRestOfParameters() const {
    std::shared_ptr<AtomFactory> atomFactory = std::make_shared<ArgonFactory>();
    std::shared_ptr<Progressbar> progressbar = std::make_shared<Progressbar>();
    this->simulation->setAtomFactory(atomFactory);
    this->simulation->setProgressbar(progressbar);
    this->simulation->setCutoffRadius(std::pow(2.0 , (1.0)/(6.0)));
    this->simulation->setRegion(this->simulation->get<Eigen::VectorXd>("initUcell") * (1.0 / std::sqrt(this->simulation->get<double>("density"))));
    this->simulation->setNumberOfAtoms(static_cast<int>(this->simulation->get<Eigen::VectorXd>("initUcell").prod()));
    this->simulation->setVelocityMagnitude(std::sqrt(this->simulation->get<int>("numberOfDimensions") *  (1.0 - (1.0/this->simulation->get<int>("numberOfAtoms"))*(this->simulation->get<double>("temperature")))));
}

void Simulation2DBuilder::initializeVectors(){
    if (!this->simulation) {
        std::cerr << "Error: simulation is not initialized!" << std::endl;
        return;
    }
    int numberOfAtoms = this->simulation->get<int>("numberOfAtoms");
    int sizeHistVel = this->simulation->get<int>("sizeHistVel");
    std::shared_ptr<AtomFactory> atomFactory = this->simulation->getAtomFactory();
    if (!atomFactory) {
        std::cerr << "Error: AtomFactory is not initialized!" << std::endl;
        return;
    }

    std::cout << "numberOfAtoms: " << numberOfAtoms << std::endl;
    std::cout << "sizeHistVel: " << sizeHistVel << std::endl;

    this->simulation->getAtoms().clear();
    this->simulation->getHistogramVelocities().clear();


    for (int i = 0; i < numberOfAtoms; ++i) {
        this->simulation->getAtoms().emplace_back(atomFactory->createAtom());
    }
    for (int i = 0; i < sizeHistVel; ++i) {
        this->simulation->getHistogramVelocities().emplace_back(0);
    }

    std::cout << "Atoms vector size after initialization: " << this->simulation->getAtoms().size() << std::endl;
    std::cout << "HistogramVelocities vector size after initialization: " << this->simulation->getHistogramVelocities().size() << std::endl;
}

void Simulation2DBuilder::initializeCoordinates() {
    auto coordinates = Eigen::VectorXd(2);
    coordinates << 0,0;
    auto gap = Eigen::VectorXd(2);
    gap << 0,0;
    //set the coordinates to two dimensions

    std::vector<std::shared_ptr<Atom2D>>& atoms = this->simulation->getAtoms();


    int nx,ny;
    // Calculate the gap as the element-wise division of region and initUcell
    gap = this->simulation->get<Eigen::VectorXd>("region").cwiseQuotient(this->simulation->get<Eigen::VectorXd>("initUcell"));

    // Reset the counter for the number of atoms initialized
    this->setNumberOfAtomIterations(0);
    for (nx = 0; nx < this->simulation->get<int>("initUcellX"); nx++) {
        for (ny = 0; ny < this->simulation->get<int>("initUcellY"); ny++) {
            // Set the coordinates to the center of the current unit cell

            coordinates(0) = (nx + 0.5);
            coordinates(1) = (ny + 0.5);

            // Scale the coordinates by the gap
            coordinates = coordinates.cwiseProduct(gap);

            // Adjust the coordinates to the region
            coordinates -= this->simulation->get<Eigen::VectorXd>("region") * 0.5;

            // Set the atom's coordinates
            atoms[this->simulation->get<int>("numberOfAtomIterations")]->setCoordinates(coordinates);

            // Increment the counter for the number of atoms initialized
            this->setNumberOfAtomIterations(this->simulation->get<int>("numberOfAtomIterations") + 1);
        }
    }
    std::cout<<"Coordinates initialized"<<std::endl;

}

void Simulation2DBuilder::initializeVelocities() {
    Eigen::VectorXd vSum = Eigen::VectorXd(2);
    Eigen::VectorXd tempVelocities = Eigen::VectorXd(2);
    vSum << 0, 0;
    tempVelocities << 0, 0;
    this->simulation->setVelocitiesSum(vSum);

    std::vector<std::shared_ptr<Atom2D>>& atoms = this->simulation->getAtoms();
    double velMag = this->simulation->get<double>("velocityMagnitude");
    int numAtoms = this->simulation->get<int>("numberOfAtoms");

    // Initialize the random number generator
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Initialize the velocities and calculate the sum of velocities
    for (auto& atom : atoms) {
        // Generate random velocity for the atom
        tempVelocities = velocity_rand(rng,dist);

        // Scale the velocity by the magnitude
        tempVelocities *= velMag;

        // Set the atom's velocity
        atom->setVelocities(tempVelocities);

        // Add the atom's velocity to the sum of velocities
        vSum += tempVelocities;
    }

    // Subtract the mean velocity to ensure the total momentum is zero
    for (auto& atom : atoms) {
        tempVelocities = atom->getVelocities();
        tempVelocities -= vSum / numAtoms;
        atom->setVelocities(tempVelocities);
    }
    std::cout<<"Velocities initialized"<<std::endl;
}

void Simulation2DBuilder::initializeAccelerations() {
    auto atoms = simulation->getAtoms();
    for (const auto& atom : atoms) {
        // Set acceleration to zero
        auto temp_accelerations = Eigen::VectorXd(2);
        temp_accelerations << 0, 0;

        atom->setAccelerations(temp_accelerations);
    }
}

void Simulation2DBuilder::setupStaticSimulation() {
    this->simulation->setCountStep(0);
    this->simulation->setCountVelocities(0);
    initializeCoordinates();
    initializeVelocities();
    initializeAccelerations();
    this->simulation->accumulateProperties(0);
}

std::shared_ptr<Simulation> Simulation2DBuilder::getSimulation() {
    return this->simulation;
}

