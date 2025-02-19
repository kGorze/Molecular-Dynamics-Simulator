//WRITTEN BY: github.com/kgorze

#ifndef LEAPFROG_H
#define LEAPFROG_H

#include <Eigen/Dense>

void leapfrog_coordinates(
    Eigen::VectorXd&            coordinates,
    double                      scale,
    const Eigen::VectorXd&      velocity
    );

void leapfrog_velocity(
    Eigen::VectorXd&        velocity,
    double                  scale,
    const Eigen::VectorXd&  acceleration
    );

#endif //LEAPFROG_H
