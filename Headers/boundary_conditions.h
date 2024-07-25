//
// Created by konrad_guest on 20/05/2024.
//

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H


#include <Eigen/Dense>

// Function to wrap a single component of a vector
inline void wrapComponent(double& component, double regionComponent) {
    if (component >= 0.5 * regionComponent) {
        component -= regionComponent;
    } else if (component < -0.5 * regionComponent) {
        component += regionComponent;
    }
}

// Function to wrap all components of a 2D vector
inline void wrapAll(Eigen::VectorXd& vec, const Eigen::VectorXd& region) {
    wrapComponent(vec.x(), region.x());
    wrapComponent(vec.y(), region.y());
}

#endif //BOUNDARY_CONDITIONS_H
