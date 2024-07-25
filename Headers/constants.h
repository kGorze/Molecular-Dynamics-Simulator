//
// Created by konrad_guest on 25/07/2024.
//

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Core>
// {x,y,z}
static const Eigen::Vector3i OFFSET_VALS[] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1},
    {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1},
    {-1,-1,1}, {0,-1,1}, {1,-1,1}
};


#endif //CONSTANTS_H
