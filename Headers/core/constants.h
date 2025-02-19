//WRITTEN BY: github.com/kgorze

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
