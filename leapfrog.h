#pragma once
#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "types.h"

void leapfrog_coordinates(VecR& coordinates, double scale, VecR& velocity ){
    coordinates.x +=  scale * velocity.x;
    coordinates.y +=  scale * velocity.y;
}

void leapfrog_velocity(VecR& velocity, double scale, VecR& accelaration){
    velocity.x +=  scale * accelaration.x;
    velocity.y +=  scale * accelaration.y;
}

#endif