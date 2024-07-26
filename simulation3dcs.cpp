//
// Created by konrad_guest on 25/07/2024.
//


#include "simulation3dcs.h"


double simulation3dcs::singleSimulationStep() {
    std::vector<std::shared_ptr<Atom>>& atoms = this->getAtoms();
    int stepCount = get<int>("stepCount");
    this->setStepCount(stepCount + 1);
    stepCount = get<int>("stepCount");
    this->setTimeNow(stepCount * get<double>("deltaT"));

    std::vector<std::tuple<int, Eigen::VectorXd>>& dataCoordinates = this->getDataCoordinates();

    clock_t begin = clock();

    leapfrogStep(1);
    applyBoundaryConditions();
    computeForces();
    leapfrogStep(2);

    evaluateProperties();
    accumulateProperties(1);

    // for (int i = 0; i < get<int>("numberOfAtoms"); ++i) {
    //     Eigen::VectorXd coordinates(3);
    //     coordinates << atoms[i]->getCoordinates().x(), atoms[i]->getCoordinates().y(), atoms[i]->getCoordinates().z();
    //     dataCoordinates.emplace_back(stepCount, coordinates);
    // }

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

void simulation3dcs::leapfrogStep(unsigned int part) {
    double deltaT = this->get<double>("deltaT");
    std::vector<std::shared_ptr<Atom>>& atoms = this->getAtoms();

    if (part == 1) {
        for (auto& atom : atoms) {
            Eigen::VectorXd velocity = atom->getVelocities();
            Eigen::VectorXd coordinates = atom->getCoordinates();
            const Eigen::VectorXd acceleration = atom->getAccelerations();

            leapfrog_velocity(velocity, 0.5 * deltaT, acceleration);
            leapfrog_coordinates(coordinates, deltaT, velocity);

            atom->setVelocities(velocity);
            atom->setCoordinates(coordinates);
            atom->setAccelerations(acceleration);
        }
    } else {
        for (auto& atom : atoms) {
            Eigen::VectorXd velocity = atom->getVelocities();
            const Eigen::VectorXd acceleration = atom->getAccelerations();

            leapfrog_velocity(velocity, 0.5 * deltaT, acceleration);

            atom->setVelocities(velocity);
            atom->setAccelerations(acceleration);
        }
    }
}

void simulation3dcs::applyBoundaryConditions() {
    Eigen::Vector2d region = this->get<Eigen::Vector2d>("region");
    auto atoms = this->getAtoms();
    for (auto& atom : atoms) {
        Eigen::Vector2d coordinates = atom->getCoordinates();
        wrap_all_Xd(coordinates, region);
        atom->setCoordinates(coordinates);
    }
}

void simulation3dcs::computeForces() {
    // Eigen::Vector3d distanceVector, inverseCellWidth, scaledPosition, periodicShift;
    // Eigen::Vector3i cellCoords, cell1, cell2;
    // std::vector<Eigen::Vector3i> neighborOffsets(std::begin(OFFSET_VALS), std::end(OFFSET_VALS));
    //
    // double forceScalar, distanceSquared, cutoffDistanceSquared, invDistanceSquared, invDistancePow6, potentialEnergy;
    // int cellIndex, molecule1, molecule2, x, y, z, neighborIndex, offset;
    // int atom1, atom2; // Declare atom1 and atom2
    //
    // cutoffDistanceSquared = std::pow(this->get<double>("cutoffRadius"), 2);
    //
    // inverseCellWidth = cells.cwiseQuotient(region);
    //
    // int totalCells = cells.x() * cells.y() * cells.z();
    // std::vector<int> cellList(numberOfAtoms + totalCells, -1);
    //
    // // Assign atoms to cells
    // for (int i = 0; i < numberOfAtoms; i++) {
    //     scaledPosition = atoms[i]->position.cwiseProduct(inverseCellWidth);
    //     cellCoords = scaledPosition.cast<int>();
    //     cellIndex = cellCoords.x() + cells.x() * (cellCoords.y() + cells.y() * cellCoords.z()) + numberOfAtoms;
    //     cellList[i] = cellList[cellIndex];
    //     cellList[cellIndex] = i;
    // }
    //
    // // Initialize accelerations
    // for (auto& atom : atoms) {
    //     atom->acceleration.setZero();
    // }
    //
    // potentialEnergySum = 0.0;
    // viralEnergySum = 0.0;
    //
    // for (z = 0; z < cells.z(); z++) {
    //     for (y = 0; y < cells.y(); y++) {
    //         for (x = 0; x < cells.x(); x++) {
    //             cell1 = Eigen::Vector3i(x, y, z);
    //             int cell1Index = cell1.x() + cells.x() * (cell1.y() + cells.y() * cell1.z()) + numberOfAtoms;
    //
    //             for (const auto& offset : neighborOffsets) {
    //                 cell2 = cell1 + offset;
    //                 periodicShift.setZero();
    //
    //                 // Apply periodic boundary conditions
    //                 for (int d = 0; d < 3; d++) {
    //                     if (cell2[d] < 0) {
    //                         cell2[d] += cells[d];
    //                         periodicShift[d] = -region[d];
    //                     } else if (cell2[d] >= cells[d]) {
    //                         cell2[d] -= cells[d];
    //                         periodicShift[d] = region[d];
    //                     }
    //                 }
    //
    //                 int cell2Index = cell2.x() + cells.x() * (cell2.y() + cells.y() * cell2.z()) + numberOfAtoms;
    //
    //                 for (atom1 = cellList[cell1Index]; atom1 >= 0; atom1 = cellList[atom1]) {
    //                     for (atom2 = cellList[cell2Index]; atom2 >= 0; atom2 = cellList[atom2]) {
    //                         if (cell1Index != cell2Index || atom2 < atom1) {
    //                             distanceVector = atoms[atom1]->position - atoms[atom2]->position + periodicShift;
    //                             distanceSquared = distanceVector.squaredNorm();
    //
    //                             if (distanceSquared < cutoffDistanceSquared) {
    //                                 invDistanceSquared = 1.0 / distanceSquared;
    //                                 invDistancePow6 = invDistanceSquared * invDistanceSquared * invDistanceSquared; // Declare and calculate invDistancePow6
    //                                 forceScalar = 48.0 * invDistancePow6 * (invDistancePow6 - 0.5) * invDistanceSquared;
    //                                 potentialEnergy = 4.0 * invDistancePow6 * (invDistancePow6 - 1.0) + 1.0;
    //
    //                                 atoms[atom1]->acceleration += forceScalar * distanceVector;
    //                                 atoms[atom2]->acceleration -= forceScalar * distanceVector;
    //
    //                                 potentialEnergySum += potentialEnergy;
    //                                 viralEnergySum += forceScalar * distanceSquared;
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
}

void simulation3dcs::evaluateProperties()
{
    std::vector<std::shared_ptr<Atom>>& atoms     = this->getAtoms();
    int numberOfAtoms                               = static_cast<int>(this->getAtoms().size());
    double velocitySquared,velocitySquaredSum       = 0; // Variable to hold the velocity squared
    Eigen::Vector2d velocitiesSum(0, 0); // Initialize the sum of velocities to zero

    for (const auto& atom : atoms)
    {
        velocitiesSum           += atom->getVelocities(); // Accumulate the sum of velocities
        velocitySquared         = atom->getVelocities().squaredNorm(); // Calculate the squared velocity
        velocitySquaredSum      += velocitySquared; // Accumulate the sum of squared velocities
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

void simulation3dcs::evaluateVelocityDistribution()
{
    int j;
    std::vector<double> hist = getHistogramVelocities();
    std::vector<std::shared_ptr<Atom>> atoms = this->getAtoms();

    if(get<int>("countVelocities") == 0){std::fill(hist.begin(), hist.end(), 0);}

    try
        {
        double deltaVelocity = static_cast<double>(get<double>("rangeVel")) / get<int>("sizeHistVel");
        if (deltaVelocity > std::numeric_limits<double>::max() || deltaVelocity < std::numeric_limits<double>::lowest()){throw mds::eeInvalidValue;}
        }
    catch(mds::ErrorCode& e)
        {
        std::cerr << "Caught exception: " << mds::getErrorCodeString(e) << '\n';
        }
    catch(...)
        {
        std::cerr << "Caught an unknown exception\n";
        }

    for (const auto& atom : atoms)
    {
        if (std::abs(deltaVelocity) > std::numeric_limits<double>::epsilon())
            {
            auto atomValue              = atom->getVelocities().norm();
            j                                 = static_cast<int>(atomValue / (std::abs(deltaVelocity)));
            hist.at(std::min(j, get<int>("sizeHistVel") - 1)) += 1;
            }
    }

    // Increment the velocity count
    setCountVelocities(get<int>("countVelocities") + 1);

    // If the count reaches the limit, normalize and print the histogram
    if(get<int>("countVelocities") == get<int>("limitVel")) {
        double histSum = std::accumulate(hist.begin(), hist.end(), 0.0);
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

void simulation3dcs::printVelocityDestribution() const
{
    std::cout << "Printing velocity distribution" << std::endl;
    for(int n = 0; n < get<int>("sizeHistVel"); n++) {

        std::cout << "Velocity: " << (n+0.5) * get<double>("rangeVel")/get<int>("sizeHistVel") << " Distribution: " <<
        histogramVelocities[n] << std::endl;
    }
};

void simulation3dcs::printSummary() const
{
    std::cout<<"\nPrinting summary"<<std::endl;
    std::cout << "Length of iteration data: " << iterationData.size() << std::endl;
    if (!iterationData.empty()) {
        const auto& data = iterationData.back(); // Get the last element
        std::cout << "Step Count: "                     << std::get<0>(data) << "\n"
                  << "Time Now: "                       << std::get<1>(data) << "\n"
                  << "Velocity Sum Value X: "           << std::get<2>(data).x() << "\n"
                  << "Velocity Sum Value Y: "           << std::get<2>(data).y() << "\n"
                  << "Potential Energy Sum: "           << std::get<3>(data) << "\n"
                  << "Potential Energy Sum Squared: "   << std::get<4>(data) << "\n"
                  << "Kinetic Energy Sum: "             << std::get<5>(data) << "\n"
                  << "Kinetic Energy Sum Squared: "     << std::get<6>(data) << "\n"
                  << "Pressure Sum: "                   << std::get<7>(data) << "\n"
                  << "Pressure Sum Squared: "           << std::get<8>(data) << "\n";
    }
}

void simulation3dcs::setIterationProperties()
{
    Eigen::Vector2d velocitySumValue    = get<Eigen::Vector2d>("velocitiesSum") / get<int>("numberOfAtoms");


    double potentialEnergySum           = get<Eigen::Vector3d>("potentialEnergy").y();
    double potentialEnergySumSquared    = get<Eigen::Vector3d>("potentialEnergy").z();
    double kineticEnergySum             = get<Eigen::Vector3d>("kineticEnergy").y();
    double kineticEnergySumSquared      = get<Eigen::Vector3d>("kineticEnergy").z();
    double pressureSum                  = get<Eigen::Vector3d>("pressure").y();
    double pressureSumSquared           = get<Eigen::Vector3d>("pressure").z();
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

void simulation3dcs::run(
    unsigned int    option
    )
{
    std::cout<<"\nStarting simulation2D"<<std::endl;
    uint64_t total_steps = get<int>("stepLimit");
    getProgressbar()->totalTicks(total_steps);

    for(int i = 0; i< total_steps; ++i){
        singleSimulationStep();
        if(i % 10 == 0) {setIterationProperties();}
        getProgressbar()->tick();
    }
    std::cout<<"\nEnding simulation2D"<<std::endl;
    //if (option == 1){
        //printSummary();
        //this->saveSummaryToFile();
        //this->saveHistogramToFile();
    //}
}

void simulation3dcs::printConfig()
{
    for (const auto& [key, value] : config) {
        std::cout << key << ": " << value << '\n';
    }
}

void simulation3dcs::plotDataHistogramVelocities(
    const std::vector<std::vector<double>>&      dataHistogramVelocities
    )
{
    auto figure = matplot::figure();
    matplot::area(dataHistogramVelocities);
    matplot::show();
}

void simulation3dcs::accumulateProperties(unsigned int code) {
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