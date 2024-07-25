//
// Created by konrad_guest on 16/05/2024.
//

#include "Headers/abstract_factory.h"

Eigen::Vector2d Atom::getCoordinates() const {
    return coordinates;
}

Eigen::Vector2d Atom::getVelocities() const {
    return velocities;
}

Eigen::Vector2d Atom::getAccelerations() const {
    return accelerations;
}

void Atom::setCoordinates(double x, double y) {
    coordinates << x, y;
}

void Atom::setVelocities(double x, double y) {
    velocities << x, y;
}

void Atom::setAccelerations(double x, double y) {
    accelerations << x, y;
}

std::shared_ptr<Atom> ArgonFactory::createAtom() {
    return std::make_shared<Argon2D>();
}
