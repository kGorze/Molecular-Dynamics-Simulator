//
// Created by konrad_guest on 16/05/2024.
//

#include "Headers/builder.h"
#include "Headers/leapfrog.h"
#include "Headers/boundary_conditions.h"
#include <iostream>
#include <random>
#include <matplot/matplot.h>

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
/*
Simulation2D::Simulation2D() {
    this->Reset();
}

void Simulation2D::Reset() {
    this->simulation = std::make_shared<Simulation>();
}
*/

void Simulation2DBuilder::ensureDirectoryExists(const std::string &path) {
    std::filesystem::path contentRoot = std::filesystem::current_path();
    std::filesystem::path dirPath = std::filesystem::path(path).parent_path();
    if (!std::filesystem::exists(dirPath)) {
        std::filesystem::create_directories(dirPath);
    }
}
template<> //PROPERTIES TEMPALTE
void Simulation2DBuilder::saveDataToFile(
    const std::string& filename,
    const std::vector<std::tuple<int, double,Eigen::Vector2d,double, double, double, double, double, double>>& data,
    unsigned int mode
    )
{

    //ensureDirectoryExists(filePath);

    //OPEN THE FILE PROTOCOL
    std::filesystem::path contentRoot = std::filesystem::current_path();
    std::filesystem::path cleanFilePath = contentRoot.remove_filename();
    std::filesystem::path relativePath = "Resources\\properties.csv";
    std::filesystem::path filePath= cleanFilePath / relativePath;


    auto& file = this->getPropertiesDataOutput();
    if(mode == 0) {
        file.open(filePath, std::ios::out | std::ios::app);
    }else if(mode == 1) {
        file.open(filePath, std::ios::out | std::ios::trunc);
    }
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing." << "\n";
    }

    // Write the header row
    file << "Step Count,Time Now,Velocity Sum Value X,Velocity Sum Value Y,Potential Energy Sum,"
         << "Potential Energy Sum Squared,Kinetic Energy Sum,Kinetic Energy Sum Squared,Pressure Sum,Pressure Sum Squared\n";

    // Write the data rows
    for (const auto& entry : data) {
        file << std::setw(5) << std::get<0>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<1>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<2>(entry).x() << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<2>(entry).y() << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<3>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<4>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<5>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<6>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<7>(entry) << ","
             << std::setw(8) << std::fixed << std::setprecision(4) << std::get<8>(entry) << "\n";
    }

    if(this->getPropertiesDataOutput().is_open()) {
        this->getPropertiesDataOutput().close();
    }

};


template <> //template for the velocities
void Simulation2DBuilder::saveDataToFile(
    const std::string& filename,
    const std::vector<std::vector<double>>& data,
    unsigned int mode)
{
    //ensureDirectoryExists(filePath);

    //OPEN THE FILE PROTOCOL
    std::filesystem::path contentRoot = std::filesystem::current_path();
    std::filesystem::path cleanFilePath = contentRoot.remove_filename();
    std::filesystem::path relativePath = "Resources\\velocity_distribution.csv";
    std::filesystem::path filePath= cleanFilePath / relativePath;

    auto& file = this->getVelocityDistributionDataOutput();
    if(mode == 0) {
        file.open(filePath, std::ios::out | std::ios::app);
    }else if(mode == 1) {
        file.open(filePath, std::ios::out | std::ios::trunc);
    }
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing." << "\n";
    }

    // Write the header row
    file << "Velocity distribution attitude\n";

    // Write the data rows
    for (const auto& row : data) {
        for (const auto& column : row) {
            file << std::setw(8) << std::fixed << std::setprecision(4) << column << ",";
        }
        file << "\n";
        file << "---\n";
    }


    if(this->getVelocityDistributionDataOutput().is_open()) {
        this->getVelocityDistributionDataOutput().close();
    }
}

template <> //COORDINATES TEMPLATE
void Simulation2DBuilder::saveDataToFile(
    const std::string& filename,
    const std::vector<std::tuple<int, double, double>>& data,
    unsigned int mode
)
{
    //ensureDirectoryExists(filePath);

    std::filesystem::path contentRoot = std::filesystem::current_path();
    std::filesystem::path cleanFilePath = contentRoot.remove_filename();
    std::filesystem::path relativePath = "Resources\\results.csv";
    std::filesystem::path filePath= cleanFilePath / relativePath;


    auto& file = this->getcoordinatesDataOutput();
    if(mode == 0) {
        file.open(filePath, std::ios::out | std::ios::app);
    }else if(mode == 1) {
        file.open(filePath, std::ios::out | std::ios::trunc);
    }
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing." << "\n";
    }

    for (const auto& entry : data) {
        file << std::setw(5) << std::get<0>(entry) << "\t"
                << std::setw(8) << std::fixed << std::setprecision(4) << std::get<1>(entry) << "\t"
                << std::setw(8) << std::fixed << std::setprecision(4) << std::get<2>(entry) << "\n";
    }

    if(this->getcoordinatesDataOutput().is_open()) {
        this->getcoordinatesDataOutput().close();
    }
}


void Simulation2DBuilder::setRestOfParameters() const {
    std::shared_ptr<AtomFactory> atomFactory = std::make_shared<ArgonFactory>();
    std::shared_ptr<Progressbar> progressbar = std::make_shared<Progressbar>();
    this->simulation->setAtomFactory(atomFactory);
    this->simulation->setProgressbar(progressbar);
    this->simulation->setCutoffRadius(std::pow(2.0 , (1.0)/(6.0)));
    this->simulation->setRegion(this->simulation->get<Eigen::Vector2d>("initUcell") * (1.0 / std::sqrt(this->simulation->get<double>("density"))));
    this->simulation->setNumberOfAtoms(static_cast<int>(this->simulation->get<Eigen::Vector2d>("initUcell").prod()));
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
    Eigen::Vector2d coordinates, gap;
    std::vector<std::shared_ptr<Atom2D>>& atoms = this->simulation->getAtoms();

    // Initialize gap and coordinates to zero
    gap(0) = 0, gap(1) =  0;
    coordinates(0) = 0; coordinates(1) = 0;

    int nx,ny;
    // Calculate the gap as the element-wise division of region and initUcell
    gap = this->simulation->get<Eigen::Vector2d>("region").cwiseQuotient(this->simulation->get<Eigen::Vector2d>("initUcell"));

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
            coordinates -= this->simulation->get<Eigen::Vector2d>("region") * 0.5;

            // Set the atom's coordinates
            atoms[this->simulation->get<int>("numberOfAtomIterations")]->setCoordinates(coordinates(0),coordinates(1));

            // Increment the counter for the number of atoms initialized
            this->setNumberOfAtomIterations(this->simulation->get<int>("numberOfAtomIterations") + 1);
        }
    }
    std::cout<<"Coordinates initialized"<<std::endl;

}

void Simulation2DBuilder::initializeVelocities() {
    Eigen::Vector2d vSum(0, 0);
    Eigen::Vector2d tempVelocities(0,0);
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
        atom->setVelocities(tempVelocities.x(), tempVelocities.y());

        // Add the atom's velocity to the sum of velocities
        vSum += tempVelocities;
    }

    // Subtract the mean velocity to ensure the total momentum is zero
    for (auto& atom : atoms) {
        tempVelocities = atom->getVelocities();
        tempVelocities -= vSum / numAtoms;
        atom->setVelocities(tempVelocities.x(), tempVelocities.y());
    }
    std::cout<<"Velocities initialized"<<std::endl;
}

void Simulation2DBuilder::initializeAccelerations() {
    auto atoms = simulation->getAtoms();
    for (const auto& atom : atoms) {
        // Set acceleration to zero
        atom->setAccelerations(0.0, 0.0);
    }
}

void Simulation2D::accumulateProperties(unsigned int code) {
    if(code == 0) {
        Eigen::Vector3d temp(0,0,0);
        this->setKineticEnergy(temp);
        this->setPotentialEnergy(temp);
        this->setPressure(temp);
    }
    else if(code == 1) {
        Eigen::Vector3d kin = this->get<Eigen::Vector3d>("kineticEnergy");
        double &kvalue = kin(0);
        this->setKineticEnergy(Eigen::Vector3d(kvalue, kvalue, std::pow(kvalue, 2)));

        Eigen::Vector3d pot = this->get<Eigen::Vector3d>("potentialEnergy");
        double &pvalue = pot(0);
        this->setPotentialEnergy(Eigen::Vector3d(pvalue, pvalue, std::pow(pvalue, 2)));

        Eigen::Vector3d pre = this->get<Eigen::Vector3d>("pressure");
        double &prvalue = pre(0);
        this->setPressure(Eigen::Vector3d(prvalue, prvalue, std::pow(prvalue, 2)));
    } else if(code == 2) {
        int stepAvg = this->get<int>("stepAvg");
        Eigen::Vector3d kin = this->get<Eigen::Vector3d>("kineticEnergy");
        double &kvalue = kin(0);
        double &ksumSquared = kin(2);
        this->setKineticEnergy(Eigen::Vector3d(kvalue, kvalue/stepAvg, std::sqrt(std::max(ksumSquared / stepAvg - std::pow(kvalue / stepAvg, 2), 0.0))));

        Eigen::Vector3d pot = this->get<Eigen::Vector3d>("potentialEnergy");
        double &pvalue = pot(0);
        double &psumSquared = pot(2);
        this->setPotentialEnergy(Eigen::Vector3d(pvalue, pvalue/stepAvg, std::sqrt(std::max(psumSquared / stepAvg - std::pow(pvalue / stepAvg, 2), 0.0))));

        Eigen::Vector3d pre = this->get<Eigen::Vector3d>("pressure");
        double &prvalue = pre(0);
        double &prsumSquared = pre(2);
        this->setPressure(Eigen::Vector3d(prvalue, prvalue/stepAvg, std::sqrt(std::max(prsumSquared / stepAvg - std::pow(prvalue / stepAvg, 2), 0.0))));
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

double Simulation2D::singleSimulationStep() {
    std::vector<std::shared_ptr<Atom2D>>& atoms = this->getAtoms();
    int stepCount = get<int>("stepCount");
    this->setStepCount(stepCount + 1);
    stepCount = get<int>("stepCount");
    this->setTimeNow(stepCount * get<double>("deltaT"));

    std::vector<std::tuple<int, double, double>>& dataCoordinates = this->getDataCoordinates();

    clock_t begin = clock();

    leapfrogStep(1);
    applyBoundaryConditions();
    computeForcesPairs();
    leapfrogStep(2);

    evaluateProperties();
    accumulateProperties(1);

    for (int i = 0; i < get<int>("numberOfAtoms"); ++i) {
        dataCoordinates.emplace_back(stepCount, atoms[i]->getCoordinates().x(), atoms[i]->getCoordinates().y());
    }

    if (stepCount % get<int>("stepAvg") == 0) {
        accumulateProperties(2);
        // if (stepCount >= get<int>("stepEquil") && (stepCount - get<int>("stepEquil")) % get<int>("stepVel") == 0) {
        // }
    }
    if(get<int>("stepCount") % 25 == 0) {
        evaluateVelocityDistribution();
        //set the histogram data to first index of dataHistogramVelocities
        std::vector<double> hist = getHistogramVelocities();
        std::vector<std::vector<double>> dataHistogramVelocities = getDataHistogramVelocities();
        dataHistogramVelocities.insert(dataHistogramVelocities.begin(), hist);
        setDataHistogramVelocities(dataHistogramVelocities);
    }
    // if (get<int>("stepCount") == get<int>("stepLimit")) {
    //     printVelocityDestribution();
    // }
    clock_t end = clock();

    double elapsed_secs = static_cast<double>(end - begin) / static_cast<double>(CLOCKS_PER_SEC);
    return elapsed_secs;
}

void Simulation2D::leapfrogStep(unsigned int part) {
    double deltaT = this->get<double>("deltaT");
    std::vector<std::shared_ptr<Atom2D>>& atoms = this->getAtoms();

    if (part == 1) {
        for (auto& atom : atoms) {
            Eigen::Vector2d velocity = atom->getVelocities();
            Eigen::Vector2d coordinates = atom->getCoordinates();
            const Eigen::Vector2d acceleration = atom->getAccelerations();

            leapfrog_velocity(velocity, 0.5 * deltaT, acceleration);
            leapfrog_coordinates(coordinates, deltaT, velocity);

            atom->setVelocities(velocity.x(), velocity.y());
            atom->setCoordinates(coordinates.x(), coordinates.y());
            atom->setAccelerations(acceleration.x(), acceleration.y());
        }
    } else {
        for (auto& atom : atoms) {
            Eigen::Vector2d velocity = atom->getVelocities();
            const Eigen::Vector2d acceleration = atom->getAccelerations();

            leapfrog_velocity(velocity, 0.5 * deltaT, acceleration);

            atom->setVelocities(velocity.x(), velocity.y());
            atom->setAccelerations(acceleration.x(), acceleration.y());
        }
    }
}

void Simulation2D::applyBoundaryConditions() {
    Eigen::Vector2d region = this->get<Eigen::Vector2d>("region");
    auto atoms = this->getAtoms();
    for (auto& atom : atoms) {
        Eigen::Vector2d coordinates = atom->getCoordinates();
        wrapAll(coordinates, region);
        atom->setCoordinates(coordinates.x(), coordinates.y());
    }
}

void Simulation2D::computeForcesPairs() {
    Eigen::Vector2d distance(0,0);
    double radiusSquared, cutoffRadiusSquared;

    auto& atoms = this->getAtoms();
    size_t numAtoms = atoms.size();
    Eigen::Vector2d region = this->get<Eigen::Vector2d>("region");
    double cutOffRadius = this->get<double>("cutOffRadius");

    cutoffRadiusSquared = cutOffRadius * cutOffRadius;

    for (auto& atom : atoms) {
        atom->setAccelerations(0,0);
    }

    this->setpotentialEnergySum(0.0);
    this->setviralEnergySum(0.0);

    for (size_t j1 = 0; j1 < numAtoms - 1; ++j1) {
        for (size_t j2 = j1 + 1; j2 < numAtoms; ++j2) {
            distance = atoms[j1]->getCoordinates() - atoms[j2]->getCoordinates();

            if (distance.x() >= 0.5 * region.x()) distance.x() -= region.x();
            else if (distance.x() < -0.5 * region.x()) distance.x() += region.x();
            if (distance.y() >= 0.5 * region.y()) distance.y() -= region.y();
            else if (distance.y() < -0.5 * region.y()) distance.y() += region.y();

            radiusSquared = distance.squaredNorm();

            if (radiusSquared < cutoffRadiusSquared) {
                double radiusSquaredInverse = 1.0 / radiusSquared;
                double radiusCubeInverse = std::pow(radiusSquaredInverse, 3);
                double forceValue = 48.0 * radiusCubeInverse * (radiusCubeInverse - 0.5) * radiusSquaredInverse;

                auto atom1 = atoms[j1]->getAccelerations();
                atom1 += forceValue * distance;
                auto atom2 = atoms[j2]->getAccelerations();
                atom2 -= forceValue * distance;
                atoms[j1]->setAccelerations(atom1.x(), atom1.y());
                atoms[j2]->setAccelerations(atom2.x(), atom2.y());

                double temp = this->get<double>("potentialEnergySum");
                temp += 4.0 * radiusCubeInverse * (radiusCubeInverse - 1.0) + 1.0;
                this->setpotentialEnergySum(temp);

                temp = this->get<double>("viralEnergySum");
                temp += forceValue * radiusSquared;
                this->setviralEnergySum(temp);
            }
        }
    }
}

void Simulation2D::evaluateProperties() {
    std::vector<std::shared_ptr<Atom2D>>& atoms = this->getAtoms();
    int numberOfAtoms = static_cast<int>(this->getAtoms().size());
    double velocitySquared; // Variable to hold the velocity squared
    Eigen::Vector2d velocitiesSum(0, 0); // Initialize the sum of velocities to zero
    double velocitySquaredSum = 0; // Initialize the sum of squared velocities to zero

    for (const auto& atom : atoms) {
        velocitiesSum += atom->getVelocities(); // Accumulate the sum of velocities
        velocitySquared = atom->getVelocities().squaredNorm(); // Calculate the squared velocity
        velocitySquaredSum += velocitySquared; // Accumulate the sum of squared velocities
    }
    this->setVelocitiesSum(velocitiesSum);

    // Calculate the kinetic energy
    Eigen::Vector3d kin = this->get<Eigen::Vector3d>("kineticEnergy");
    kin.x() = ((0.5 * velocitySquaredSum)/numberOfAtoms );
    this->setKineticEnergy(kin);

    // Calculate the total energy
    Eigen::Vector3d pot = this->get<Eigen::Vector3d>("potentialEnergy");
    pot.x() = kin.x() + (this->get<double>("potentialEnergySum")/ numberOfAtoms );
    this->setPotentialEnergy(pot);

    // Calculate the pressure(2 for 2D simulations)
    Eigen::Vector3d pre = this->get<Eigen::Vector3d>("pressure");
    pre.x() = this->get<double>("density") * (velocitySquaredSum - this->get<double>("viralEnergySum"))/(numberOfAtoms * 2);
    this->setPressure(pre);
}

void Simulation2D::evaluateVelocityDistribution() {
    double deltaVelocity = get<double>("deltaVelocity");
    deltaVelocity = 0.0;
    int j;

    double histSum = get<double>("histogramSum");
    histSum = 0.0;


    std::vector<double> hist = getHistogramVelocities();
    std::vector<std::shared_ptr<Atom2D>> atoms = this->getAtoms();

    // If the velocity count is zero, initialize the histogram array to zeros
    if(get<int>("countVelocities") == 0) {
        std::fill(hist.begin(), hist.end(), 0);
    }

    // Calculate the bin width for the histogram
    deltaVelocity = (get<double>("rangeVel") / get<int>("sizeHistVel"));

    // Determine the histogram bin index for the current molecule's velocity magnitude
    for (const auto& atom : atoms) {
        if (std::abs(deltaVelocity) > std::numeric_limits<double>::epsilon()) {
            auto atomValue = atom->getVelocities().norm();
            j = static_cast<int>(atomValue / (std::abs(deltaVelocity)));
            // Increment the histogram bin, using Min to ensure the index does not exceed the array bounds
            hist.at(std::min(j, get<int>("sizeHistVel") - 1)) += 1;
        }
    }

    // Increment the velocity count
    setCountVelocities(get<int>("countVelocities") + 1);

    // If the count reaches the limit, normalize and print the histogram
    if(get<int>("countVelocities") == get<int>("limitVel")) {
        setHistogramSum(std::accumulate(hist.begin(), hist.end(), 0.0));

        // Calculate the total sum of the histogram bins
        for(auto& value: hist) {
            histSum += value;
        }

        // Normalize the histogram bins by dividing by the total sum
        for (auto& value : hist) {
            value /= histSum;
        }

        // Reset the velocity count
        setCountVelocities(0);
    }

    setHistogramVelocities(hist);

    setEntropyFunction(0.0);
    for(int i = 0; i < get<int>("sizeHistVel"); i++) {
        if(hist[i] > 0) {
            double entf = get<double>("entropyFunction");
            entf += hist[i] * std::log(hist[i]) / ((j+0.5)*get<double>("deltaVelocity"));
            setEntropyFunction(entf);
        }
    }
}

void Simulation2D::printVelocityDestribution() const {
    std::cout << "Printing velocity distribution" << std::endl;
    for(int n = 0; n < get<int>("sizeHistVel"); n++) {

        std::cout << "Velocity: " << (n+0.5) * get<double>("rangeVel")/get<int>("sizeHistVel") << " Distribution: " <<
        histogramVelocities[n] << std::endl;
    }
};

void Simulation2D::printSummary() const {
    std::cout<<"\nPrinting summary"<<std::endl;
    std::cout << "Length of iteration data: " << iterationData.size() << std::endl;
    if (!iterationData.empty()) {
        const auto& data = iterationData.back(); // Get the last element
        std::cout << "Step Count: " << std::get<0>(data) << "\n"
                  << "Time Now: " << std::get<1>(data) << "\n"
                  << "Velocity Sum Value X: " << std::get<2>(data).x() << "\n"  // Access the x component
                  << "Velocity Sum Value Y: " << std::get<2>(data).y() << "\n"  // Access the y component
                  << "Potential Energy Sum: " << std::get<3>(data) << "\n"
                  << "Potential Energy Sum Squared: " << std::get<4>(data) << "\n"
                  << "Kinetic Energy Sum: " << std::get<5>(data) << "\n"
                  << "Kinetic Energy Sum Squared: " << std::get<6>(data) << "\n"
                  << "Pressure Sum: " << std::get<7>(data) << "\n"
                  << "Pressure Sum Squared: " << std::get<8>(data) << "\n";
    }
}

void Simulation2D::setIterationProperties() {
    Eigen::Vector2d velocitySumValue = get<Eigen::Vector2d>("velocitiesSum") / get<int>("numberOfAtoms");


    double potentialEnergySum = get<Eigen::Vector3d>("potentialEnergy").y();
    double potentialEnergySumSquared = get<Eigen::Vector3d>("potentialEnergy").z();
    double kineticEnergySum = get<Eigen::Vector3d>("kineticEnergy").y();
    double kineticEnergySumSquared = get<Eigen::Vector3d>("kineticEnergy").z();
    double pressureSum = get<Eigen::Vector3d>("pressure").y();
    double pressureSumSquared = get<Eigen::Vector3d>("pressure").z();
    iterationData.emplace_back(
        get<int>("stepCount"),
        get<double>("timeNow"),
        velocitySumValue,
        potentialEnergySum,
        potentialEnergySumSquared,
        kineticEnergySum,
        kineticEnergySumSquared,
        pressureSum,
        pressureSumSquared
        );

};

void Simulation2D::run(unsigned int option) {
    //fix the option in run function
    std::cout << "Starting simulation2D" << std::endl;
    uint64_t total_steps = get<int>("stepLimit");
    getProgressbar()->totalTicks(total_steps);

    while(true) {
        singleSimulationStep();
        if(this->countStep%10 == 0) {
            setIterationProperties();
        }
        if (get<int>("stepCount") >= total_steps) {
            break;
        }
        getProgressbar()->tick();
    }
    std::cout<<"\nEnding simulation2D"<<std::endl;
    if (option == 1){
        //printSummary();
        //this->saveSummaryToFile();
        //this->saveHistogramToFile();
    }
}


/* void Simulation2D::setAtomSize(const int x_axis,const int y_axis){
    this->atoms.assign(x_axis*y_axis, Atom2D());
}

void Simulation2D::initializeAtoms() {
    if(this->atoms.empty()) {
        std::cout << "Atoms not initialized" << std::endl;
        return;
    }
    for(int i = 0; i< this->atoms.size(); i++) {
        this->atoms[i].setCoordinates(i,i);
        simulation->atoms[i].velocities = std::make_pair(0, 0);
        simulation->atoms[i].accelerations = std::make_pair(0, 0);
    }
    std::cout << "Atoms initialized" << std::endl;
}*/

/*void Simulation2D::printAtomsSize() const {
    std::cout << "Atoms size: " << atoms.size() << std::endl;
} */

/* void Simulation2D::printAtoms() const {
    for (int i = 0; i < atoms.size()/2; i++) {
        vector2d temp = atoms[i].getCoordinates();
        vector2d temp2 = atoms[i].getVelocities();
        vector2d temp3 = atoms[i].getAccelerations();
        std::cout << "Atom " << i << " coordinates: " << temp.first<< ", "<< temp.second << std::endl;
        std::cout << "Atom " << i << " velocities: " << temp2.first<< ", "<< temp2.second << std::endl;
        std::cout << "Atom " << i << " accelerations: " << temp3.first<< ", "<< temp3.second << std::endl;
    }
} */

void Simulation2D::printConfig(){
    for (const auto& pair : config) {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }
}

std::shared_ptr<Simulation> Simulation2DBuilder::getSimulation() {
    return this->simulation;
}

void Simulation2D::plotDataHistogramVelocities(const std::vector<std::vector<double>>& dataHistogramVelocities) {
    // Create a new figure
    auto figure = matplot::figure();

    // Plot the data as an area plot
    matplot::area(dataHistogramVelocities);

    // Display the figure
    matplot::show();
}

