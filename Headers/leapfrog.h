//
// Created by konrad_guest on 19/05/2024.
//

#ifndef LEAPFROG_H
#define LEAPFROG_H

#include <Eigen/Dense>

void leapfrog_coordinates(Eigen::Vector2d& coordinates, double scale, const Eigen::Vector2d& velocity) {
    coordinates += scale * velocity;
}

void leapfrog_velocity(Eigen::Vector2d& velocity, double scale, const Eigen::Vector2d& acceleration) {
    velocity += scale * acceleration;
}

#endif //LEAPFROG_H
