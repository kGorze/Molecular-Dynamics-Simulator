#ifndef RAND2D_H
#define RAND2D_H

#include <random>
#include "types.h"


void velocity_rand(Mol *mol){
    mol->velocity.x = std::cos(2 * M_PI * static_cast<double>(rand()) / RAND_MAX);
    mol->velocity.y = std::cos(2 * M_PI * static_cast<double>(rand()) / RAND_MAX);
}

#endif // RAND2D_H
