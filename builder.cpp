//
// Created by konrad_guest on 16/05/2024.
//

#include "Headers/builder.h"
#include <iostream>
#include <random>

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

    // Get current date and time
    auto now      = std::chrono::system_clock::now();
    std::time_t now_time                = std::chrono::system_clock::to_time_t(now);
    std::string date_time               = std::ctime(&now_time);
    date_time.pop_back();

    date_time = date_time.substr(0, date_time.size()-6);

    std::replace(date_time.begin(), date_time.end(), ':', '_');
    std::replace(date_time.begin(), date_time.end(), ' ', '_');
    std::replace(date_time.begin(), date_time.end(), '.', '_');


    // Get number of atoms, iterations, and density
    std::string num_atoms       = std::to_string(simulation->get<int>("numberOfAtoms"));
    std::string iterations      = std::to_string(simulation->get<int>("stepLimit"));
    std::string density         = std::to_string(simulation->get<double>("density"));

    // Get the current working directory
    std::string current_path = std::filesystem::current_path().string();


    // Construct directory name
    std::string dir_name = current_path + "\\Resources\\" + date_time + "_" + num_atoms + "_" + iterations + "_" + density;

    // Create directory
    std::filesystem::path dir_path(dir_name);
    if (!std::filesystem::exists(dir_path))
    {
        std::filesystem::create_directories(dir_path);
        std::cout << "Directory created: " << dir_path << std::endl;
    }

    //OPEN THE FILE PROTOCOL
    std::string relativePath      = "properties.csv";

    // Construct file path
    std::filesystem::path filePath = dir_name + "\\" + relativePath;


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
    const std::string&                      filename,
    const std::vector<std::vector<double>>& data,
    unsigned int                            mode
    )
{
    //ensureDirectoryExists(filePath);

    // Get current date and time
    auto now      = std::chrono::system_clock::now();
    std::time_t now_time                = std::chrono::system_clock::to_time_t(now);
    std::string date_time               = std::ctime(&now_time);
    date_time.pop_back();

    date_time = date_time.substr(0, date_time.size()-6);


    std::replace(date_time.begin(), date_time.end(), ':', '_');
    std::replace(date_time.begin(), date_time.end(), ' ', '_');
    std::replace(date_time.begin(), date_time.end(), '.', '_');


    // Get number of atoms, iterations, and density
    std::string num_atoms       = std::to_string(simulation->get<int>("numberOfAtoms"));
    std::string iterations      = std::to_string(simulation->get<int>("stepLimit"));
    std::string density         = std::to_string(simulation->get<double>("density"));

    // Get the current working directory
    std::string current_path = std::filesystem::current_path().string();


    // Construct directory name
    std::string dir_name = current_path + "\\Resources\\" + date_time + "_" + num_atoms + "_" + iterations + "_" + density;

    // Create directory
    std::filesystem::path dir_path(dir_name);
    if (!std::filesystem::exists(dir_path))
    {
        std::filesystem::create_directories(dir_path);
        std::cout << "Directory created: " << dir_path << std::endl;
    }

    std::string relativePath      = "velocity_distribution.csv";

    // Construct file path
    std::filesystem::path filePath = dir_name + "\\" + relativePath;


    auto& file = this->getVelocityDistributionDataOutput();
    if(mode == 0)
    {
        file.open(filePath, std::ios::out | std::ios::app);
    }else if(mode == 1)
    {
        file.open(filePath, std::ios::out | std::ios::trunc);
    }
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file for writing." << "\n";
    }

    // Write the header row
    file << "Velocity distribution attitude\n";

    // Write the data rows
    for (const auto& row : data) {
        for (const auto& column : row)
        {
            file << std::setw(8) << std::fixed << std::setprecision(4) << column << ",";
        }
            file << "\n"<< "---\n";
    }

    if(this->getVelocityDistributionDataOutput().is_open()) {this->getVelocityDistributionDataOutput().close();}
}

template <> //COORDINATES TEMPLATE
void Simulation2DBuilder::saveDataToFile(
    const std::string&                                  filename,
    const std::vector<std::tuple<int, double, double>>& data,
    unsigned int                                        mode
)
{
    //ensureDirectoryExists(filePath);


    auto now      = std::chrono::system_clock::now();
    std::time_t now_time                = std::chrono::system_clock::to_time_t(now);
    std::string date_time               = std::ctime(&now_time);
    date_time.pop_back();

    date_time = date_time.substr(0, date_time.size()-6);


    std::replace(date_time.begin(), date_time.end(), ':', '_');
    std::replace(date_time.begin(), date_time.end(), ' ', '_');
    std::replace(date_time.begin(), date_time.end(), '.', '_');


    // Get number of atoms, iterations, and density
    std::string num_atoms       = std::to_string(simulation->get<int>("numberOfAtoms"));
    std::string iterations      = std::to_string(simulation->get<int>("stepLimit"));
    std::string density         = std::to_string(simulation->get<double>("density"));

    // Get the current working directory
    std::string current_path = std::filesystem::current_path().string();


    // Construct directory name
    std::string dir_name = current_path + "\\Resources\\" + date_time + "_" + num_atoms + "_" + iterations + "_" + density;

    // Create directory
    std::filesystem::path dir_path(dir_name);
    if (!std::filesystem::exists(dir_path))
    {
        std::filesystem::create_directories(dir_path);
        std::cout << "Directory created: " << dir_path << std::endl;
    }
    std::string relativePath      = "results.csv";


    // Construct file path
    std::filesystem::path filePath = dir_name + "\\" + relativePath;


    auto& file = this->getcoordinatesDataOutput();
    if(mode == 0)
    {
        file.open(filePath, std::ios::out | std::ios::app);
    }else if(mode == 1)
    {
        file.open(filePath, std::ios::out | std::ios::trunc);
    }
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file for writing." << "\n";
    }

    for (const auto& entry : data) {
        file    << std::setw(5) << std::get<0>(entry) << "\t"
                << std::setw(8) << std::fixed << std::setprecision(4) << std::get<1>(entry) << "\t"
                << std::setw(8) << std::fixed << std::setprecision(4) << std::get<2>(entry) << "\n";
    }

    if(this->getcoordinatesDataOutput().is_open()) {this->getcoordinatesDataOutput().close();}
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

std::shared_ptr<Simulation> Simulation2DBuilder::getSimulation() {
    return this->simulation;
}

