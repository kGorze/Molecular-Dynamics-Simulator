//
// Created by konrad_guest on 16/05/2024.
//

#include "Headers/abstract_factory.h"

Eigen::Vector2d Atom2D::getCoordinates() const {
    return coordinates;
}

Eigen::Vector2d Atom2D::getVelocities() const {
    return velocities;
}

Eigen::Vector2d Atom2D::getAccelerations() const {
    return accelerations;
}

void Atom2D::setCoordinates(double x, double y) {
    coordinates << x, y;
}

void Atom2D::setVelocities(double x, double y) {
    velocities << x, y;
}

void Atom2D::setAccelerations(double x, double y) {
    accelerations << x, y;
}

std::shared_ptr<Atom2D> ArgonFactory::createAtom() {
    return std::make_shared<Argon2D>();
}
