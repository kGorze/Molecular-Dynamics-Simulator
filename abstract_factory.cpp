//
// Created by konrad_guest on 16/05/2024.
//

#include "Headers/abstract_factory.h"

Eigen::VectorXd Atom::getCoordinates() const {
    return coordinates;
}

Eigen::VectorXd Atom::getVelocities() const {
    return velocities;
}

Eigen::VectorXd Atom::getAccelerations() const {
    return accelerations;
}

void Atom::setCoordinates(const Eigen::VectorXd& coords) {
    coordinates << coords.x(), coords.y();
}

void Atom::setVelocities(const Eigen::VectorXd& vels) {
    velocities << vels.x(), vels.y();
}

void Atom::setAccelerations(const Eigen::VectorXd& accels) {
    accelerations << accels.x(), accels.y();
}

std::shared_ptr<Atom> ArgonFactory::createAtom() {
    return std::make_shared<Argon2D>();
}
