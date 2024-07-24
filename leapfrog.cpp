//
// Created by konrad_guest on 25/07/2024.
//
#include "leapfrog.h"

void leapfrog_coordinates(
    Eigen::Vector2d&            coordinates,
    double                      scale,
    const Eigen::Vector2d&      velocity
    )
{
    coordinates += scale * velocity;
}

void leapfrog_velocity(
    Eigen::Vector2d&        velocity,
    double                  scale,
    const Eigen::Vector2d&  acceleration
    )
{
    velocity += scale * acceleration;
}