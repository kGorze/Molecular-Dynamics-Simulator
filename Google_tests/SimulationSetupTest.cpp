//WRITTEN BY: github.com/kgorze

#include "gtest/gtest.h"
#include "Headers/simulation/simulation2dap.h"
#include "Headers/simulation/builder.h"
#include "Headers/utils/progressbar.h"


class Simulation2DAPFixture : public ::testing::Test {
protected:
    virtual void SetUp()
    {
        std::shared_ptr<AtomFactory> atomFactory = std::make_shared<ArgonFactory>();
        std::shared_ptr<Progressbar> progressbar = std::make_shared<Progressbar>();

        Eigen::VectorXd Ucell = Eigen::VectorXd(2);
        Ucell << 5,5;

        sim = new simulation2dap();

        // Initialize simulation parameters here
        sim->setDeltaT(0.005);
        sim->setDensity(0.5);
        sim->setInitUcell(Ucell);
        sim->setStepAvg(50);
        sim->setStepEquil(0);
        sim->setStepLimit(100);
        sim->setTemperature(1.0);
        sim->setLimitVel(4);
        sim->setRangeVel(3.0);
        sim->setSizeHistVel(50);
        sim->setStepVel(5);
        sim->setRandSeed(18);
        sim->setNumberOfDimensions(2);

        sim->setAtomFactory(atomFactory);
        sim->setProgressbar(progressbar);
        sim->setCutoffRadius(std::pow(2.0, (1.0) / (6.0)));

        sim->setRegion(Ucell * (1.0 / std::sqrt(0.5)));
        sim->setNumberOfAtoms(Ucell.prod());
        sim->setVelocityMagnitude(std::sqrt(2 * (1.0 - (1.0 / Ucell.prod()) * (1.0))));

        //initialzie vectors
        for(int i = 0; i< Ucell.prod(); ++i){
            sim->getAtoms().emplace_back(atomFactory->createAtom());
        }

        sim->setCountStep(0);
        sim->setCountVelocities(0);

        Eigen::VectorXd coordinates = Eigen::VectorXd(2);
        coordinates << 0,0;
        Eigen::VectorXd gap = Eigen::VectorXd(2);
        gap << 0,0;

        int nx,ny;
        std::vector<std::shared_ptr<Atom2D>>& atoms = sim->getAtoms();
        gap = sim->get<Eigen::VectorXd>("region").cwiseQuotient(sim->get<Eigen::VectorXd>("initUcell"));
        sim->setNumberOfAtomIterations(0);
        for (nx = 0; nx < sim->get<int>("initUcellX"); nx++) {
            for (ny = 0; ny < sim->get<int>("initUcellY"); ny++) {
                coordinates(0) = (nx + 0.5);
                coordinates(1) = (ny + 0.5);
                coordinates = coordinates.cwiseProduct(gap);
                coordinates -= sim->get<Eigen::VectorXd>("region") * 0.5;
                atoms[sim->get<int>("numberOfAtomIterations")]->setCoordinates(coordinates);
                sim->setNumberOfAtomIterations(sim->get<int>("numberOfAtomIterations") + 1);
            }
        }

        Eigen::VectorXd vSum = Eigen::VectorXd(2);
        Eigen::VectorXd tempVelocities = Eigen::VectorXd(2);
        vSum << 0, 0;
        tempVelocities << 0, 0;
        sim->setVelocitiesSum(vSum);

        atoms = sim->getAtoms();
        double velMag = sim->get<double>("velocityMagnitude");
        int numAtoms = sim->get<int>("numberOfAtoms");

        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (auto& atom : atoms) {
            tempVelocities = velocity_rand(rng,dist);
            tempVelocities *= velMag;
            atom->setVelocities(tempVelocities);
            vSum += tempVelocities;
        }
        for (auto& atom : atoms) {
            tempVelocities = atom->getVelocities();
            tempVelocities -= vSum / numAtoms;
            atom->setVelocities(tempVelocities);
        }

        atoms = sim->getAtoms();
        for (const auto& atom : atoms) {
            auto temp_accelerations = Eigen::VectorXd(2);
            temp_accelerations << 0, 0;
            atom->setAccelerations(temp_accelerations);
        }
        sim->accumulateProperties(0);
    }

    virtual void TearDown() {
        delete sim;
        sim = nullptr;
    }

    simulation2dap* sim;
};

TEST_F(Simulation2DAPFixture, ConfigurationSetup) {
    EXPECT_DOUBLE_EQ(sim->get<double>("deltaT"), 0.005);
    EXPECT_DOUBLE_EQ(sim->get<double>("density"), 0.5);
    EXPECT_EQ(sim->get<int>("initUcellX"), 5);
    EXPECT_EQ(sim->get<int>("initUcellY"), 5);
    EXPECT_EQ(sim->get<int>("stepAvg"), 50);
    EXPECT_EQ(sim->get<int>("stepEquil"), 0);
    EXPECT_EQ(sim->get<int>("stepLimit"), 100);
    EXPECT_DOUBLE_EQ(sim->get<double>("temperature"), 1.0);
    EXPECT_EQ(sim->get<int>("limitVel"), 4);
    EXPECT_DOUBLE_EQ(sim->get<double>("rangeVel"), 3.0);
    EXPECT_EQ(sim->get<int>("sizeHistVel"), 50);
    EXPECT_EQ(sim->get<int>("stepVel"), 5);
    EXPECT_EQ(sim->get<int>("randSeed"), 18);
    EXPECT_EQ(sim->get<int>("numberOfDimensions"), 2);
    EXPECT_DOUBLE_EQ(sim->get<double>("cutOffRadius"), std::pow(2.0, (1.0) / (6.0)));
}

TEST_F(Simulation2DAPFixture, AtomInitialization) {
    EXPECT_EQ(sim->getAtoms().size(), sim->get<int>("numberOfAtoms"));
    EXPECT_EQ(sim->get<int>("numberOfAtoms"), 25);  // 5x5 grid
}


TEST_F(Simulation2DAPFixture, SingleSimulationStep)
{
    EXPECT_EQ(sim->get<int>("stepCount"), 0);
    double elapsed_time = sim->singleSimulationStep();
    EXPECT_GT(elapsed_time, 0.0);
    EXPECT_EQ(sim->get<int>("stepCount"), 1);
}

TEST_F(Simulation2DAPFixture, LeapfrogStep)
{
    std::vector<Eigen::Vector2d> initialCoordinates;
    std::vector<Eigen::Vector2d> initialVelocities;
    for (const auto& atom : sim->getAtoms()) {
        initialCoordinates.push_back(atom->getCoordinates());
        initialVelocities.push_back(atom->getVelocities());
    }

    // Perform leapfrog step
    sim->leapfrogStep(1);

    // Check that coordinates and velocities have changed
    for (size_t i = 0; i < sim->getAtoms().size(); ++i) {
        EXPECT_NE(sim->getAtoms()[i]->getCoordinates(), initialCoordinates[i]);
        EXPECT_NE(sim->getAtoms()[i]->getVelocities(), initialVelocities[i]);
    }
}


TEST_F(Simulation2DAPFixture, ApplyBoundaryConditions)
{
    // Set up atoms outside the boundary
    sim->getAtoms()[0]->setCoordinates(Eigen::Vector2d(15.0, 15.0));

    // Apply boundary conditions
    sim->applyBoundaryConditions();

    // Check that the atom is now inside the boundary
    Eigen::Vector2d coords = sim->getAtoms()[0]->getCoordinates();
    EXPECT_LT(coords.x(), sim->get<Eigen::VectorXd>("region")(0));
    EXPECT_LT(coords.y(), sim->get<Eigen::VectorXd>("region")(1));
    EXPECT_GE(coords.x(), 0.0);
    EXPECT_GE(coords.y(), 0.0);
}

TEST_F(Simulation2DAPFixture, ComputeForces)
{
    for (auto& atom : sim->getAtoms()) {
        atom->setAccelerations(Eigen::Vector2d::Zero());
    }

    // Compute forces
    sim->computeForces();

    // Check that forces (accelerations) have been computed
    bool allZero = true;
    for (const auto& atom : sim->getAtoms()) {
        if (!atom->getAccelerations().isZero()) {
            allZero = false;
            break;
        }
    }
    EXPECT_FALSE(allZero);
}

TEST_F(Simulation2DAPFixture, EvaluateProperties) {
    sim->evaluateProperties();

    EXPECT_GT(sim->get<Eigen::Vector3d>("kineticEnergy").norm(), 0.0);
    EXPECT_GT(sim->get<Eigen::Vector3d>("potentialEnergy").norm(), 0.0);
    EXPECT_GT(sim->get<Eigen::Vector3d>("pressure").norm(), 0.0);
}

TEST_F(Simulation2DAPFixture, EvaluateVelocityDistribution) {
    sim->evaluateVelocityDistribution();

    const auto& histogram = sim->getHistogramVelocities();
    EXPECT_FALSE(histogram.empty());
    EXPECT_EQ(histogram.size(), sim->get<int>("sizeHistVel"));
}

TEST_F(Simulation2DAPFixture, AccumulateProperties) {
    int initialCount = sim->get<int>("countStep");
    sim->accumulateProperties(0);
    EXPECT_EQ(sim->get<int>("countStep"), initialCount + 1);
}