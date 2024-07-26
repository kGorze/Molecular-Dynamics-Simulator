//
// Created by konrad_guest on 20/05/2024.
//

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H


#include <Eigen/Dense>

void wrap_component_2d(double& component, double regionComponent);
void wrap_all_2d(Eigen::VectorXd& vec, const Eigen::VectorXd& region);


void wrap_component_Xd(double& component, double regionComponent);
void wrap_all_Xd(Eigen::Vector2d& vec, const Eigen::Vector2d& region);



#endif //BOUNDARY_CONDITIONS_H
