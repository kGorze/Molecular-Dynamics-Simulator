#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "types.h"

void leapfrog_coordinates(VecR* coordinates, float scale, VecR* velocity ){
    coordinates->x +=  scale * velocity->x;
    coordinates->y +=  scale * velocity->y;
}

void leapfrog_velocity(VecR* velocity, float scale, VecR* accelaration){
    velocity->x +=  scale * accelaration->x;
    velocity->y +=  scale * accelaration->y;
}

#endif // LEAPFROG_H