//WRITTEN BY: github.com/kgorze

#include "Headers/simulation/abstract_factory.h"

Eigen::VectorXd Atom2D::getCoordinates() const {
    return coordinates;
}

Eigen::VectorXd Atom2D::getVelocities() const {
    return velocities;
}

Eigen::VectorXd Atom2D::getAccelerations() const {
    return accelerations;
}

void Atom2D::setCoordinates(const Eigen::VectorXd& coords) {
    coordinates << coords.x(), coords.y();
}

void Atom2D::setVelocities(const Eigen::VectorXd& vels) {
    velocities << vels.x(), vels.y();
}

void Atom2D::setAccelerations(const Eigen::VectorXd& accels) {
    accelerations << accels.x(), accels.y();
}


std::shared_ptr<Atom2D> ArgonFactory::createAtom() {
    return std::make_shared<Atom2D>();
}
