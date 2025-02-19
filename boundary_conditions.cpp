//WRITTEN BY: github.com/kgorze

#include "Headers/simulation/boundary_conditions.h"

// Function to wrap a single component of a vector
void wrap_component_2d(double& component, double regionComponent) {
    if (component >= 0.5 * regionComponent) {
        component -= regionComponent;
    } else if (component < -0.5 * regionComponent) {
        component += regionComponent;
    }
}

// Function to wrap all components of a 2D vector
void wrap_all_2d(Eigen::VectorXd& vec, const Eigen::VectorXd& region) {
    wrap_component_2d(vec.x(), region.x());
    wrap_component_2d(vec.y(), region.y());
}

void wrap_component_Xd(double& component, double regionComponent){}
void wrap_all_Xd(Eigen::Vector2d& vec, const Eigen::Vector2d& region){}