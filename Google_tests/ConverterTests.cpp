//
// Created by konrad_guest on 27/07/2024.
//
#include "gtest/gtest.h"
#include "simulation2dap.h"


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
        // Add more parameters as needed
    }

    virtual void TearDown() {
        delete sim;
    }

    simulation2dap* sim;
};

TEST_F(Simulation2DAPFixture, SingleSimulationStep)
{
    double elapsed_time = sim->singleSimulationStep();
    EXPECT_GT(elapsed_time, 0.0);
    EXPECT_EQ(sim->get<int>("stepCount"), 1);
}

TEST_F(Simulation2DAPFixture, LeapfrogStep)
{
    // Set up initial conditions
    std::vector<std::shared_ptr<Atom2D>> atoms;
    for (int i = 0; i < 10; ++i) {
        auto atom = std::make_shared<Atom2D>();
        atom->setCoordinates(Eigen::Vector2d::Random());
        atom->setVelocities(Eigen::Vector2d::Random());
        atom->setAccelerations(Eigen::Vector2d::Zero());
        atoms.push_back(atom);
    }
    sim->setAtoms(atoms);

    // Perform leapfrog step
    sim->leapfrogStep(1);

    // Check that coordinates and velocities have changed
    for (const auto& atom : sim->getAtoms()) {
        EXPECT_FALSE(atom->getCoordinates().isZero());
        EXPECT_FALSE(atom->getVelocities().isZero());
    }
}


TEST_F(Simulation2DAPFixture, ApplyBoundaryConditions)
{
    // Set up atoms outside the boundary
    std::vector<std::shared_ptr<Atom2D>> atoms;
    auto atom = std::make_shared<Atom2D>();
    atom->setCoordinates(Eigen::Vector2d(15.0, 15.0));
    atoms.push_back(atom);
    sim->setAtoms(atoms);

    // Apply boundary conditions
    sim->applyBoundaryConditions();

    // Check that the atom is now inside the boundary
    Eigen::Vector2d coords = sim->getAtoms()[0]->getCoordinates();
    EXPECT_LT(coords.x(), 10.0);
    EXPECT_LT(coords.y(), 10.0);
    EXPECT_GE(coords.x(), 0.0);
    EXPECT_GE(coords.y(), 0.0);
}

TEST_F(Simulation2DAPFixture, ComputeForces)
{
    // Set up atoms
    std::vector<std::shared_ptr<Atom2D>> atoms;
    for (int i = 0; i < 2; ++i) {
        auto atom = std::make_shared<Atom2D>();
        atom->setCoordinates(Eigen::Vector2d::Random() * 5.0);
        atom->setAccelerations(Eigen::Vector2d::Zero());
        atoms.push_back(atom);
    }
    sim->setAtoms(atoms);

    // Compute forces
    sim->computeForces();

    // Check that forces (accelerations) have been computed
    for (const auto& atom : sim->getAtoms()) {
        EXPECT_FALSE(atom->getAccelerations().isZero());
    }
}
