//WRITTEN BY: github.com/kgorze


#include "Headers/utils/leapfrog.h"

void leapfrog_coordinates(
    Eigen::VectorXd&            coordinates,
    double                      scale,
    const Eigen::VectorXd&      velocity
    )
{
    coordinates += scale * velocity;
}

void leapfrog_velocity(
    Eigen::VectorXd&        velocity,
    double                  scale,
    const Eigen::VectorXd&  acceleration
    )
{
    velocity += scale * acceleration;
}