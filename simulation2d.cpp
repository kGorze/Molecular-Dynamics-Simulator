//
// Created by konrad_guest on 24/07/2024.
//

#include "simulation2d.h"

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
        try {
            std::vector<double> hist = getHistogramVelocities();
            if(hist.size() == 0) {
                throw mds::eeInvalidValue;
            }
            std::vector<std::vector<double>> dataHistogramVelocities = getDataHistogramVelocities();
            dataHistogramVelocities.insert(dataHistogramVelocities.begin(), hist);
            setDataHistogramVelocities(dataHistogramVelocities);
        } catch(mds::ErrorCode& e) {
            std::cerr << "Caught exception: " << mds::getErrorCodeString(e) << '\n';
        } catch(...) {
            std::cerr << "Caught an unknown exception\n";
        }

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
    try {
        double temp = static_cast<double>(get<double>("rangeVel")) / get<int>("sizeHistVel");
        if (temp > std::numeric_limits<double>::max() || temp < std::numeric_limits<double>::lowest()) {
            throw mds::eeInvalidValue;
        }
        deltaVelocity = temp;
    } catch(mds::ErrorCode& e) {
        std::cerr << "Caught exception: " << mds::getErrorCodeString(e) << '\n';
    } catch(...) {
        std::cerr << "Caught an unknown exception\n";
    }


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
        for(auto& value: hist) {histSum += value;}
        // Normalize the histogram bins by dividing by the total sum
        for (auto& value : hist) {value /= histSum;}

        // Reset the velocity count
        setCountVelocities(0);
    }

    setHistogramVelocities(hist);

    setEntropyFunction(0.0);
    for(int i = 0; i < get<int>("sizeHistVel"); i++) {
        if(hist[i] > 0) {
            double entf = get<double>("entropyFunction");

            try {
                double temp = 0.0;
                if(hist[i] > 1 and get<double>("deltaVelocity") > 0.0){
                    temp = (hist[i] * std::log(hist[i])) / ((j+0.5)*get<double>("deltaVelocity"));
                }
                else {
                    temp = hist[i] / 0.01;
                }

                if (temp > std::numeric_limits<double>::max() || temp < std::numeric_limits<double>::lowest()) {
                    throw mds::eeInvalidValue;
                }
                entf += temp;
            } catch(mds::ErrorCode& e) {
                std::cerr << "Caught exception: " << mds::getErrorCodeString(e) << '\n';
            } catch(...) {
                std::cerr << "Caught an unknown exception\n";
            }


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



void Simulation2D::plotDataHistogramVelocities(const std::vector<std::vector<double>>& dataHistogramVelocities) {
    // Create a new figure
    auto figure = matplot::figure();

    // Plot the data as an area plot
    matplot::area(dataHistogramVelocities);

    // Display the figure
    matplot::show();
}