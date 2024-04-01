#pragma once
#ifndef RAND2D_H
#define RAND2D_H

#include "types.h"

# define M_PI           3.14159265358979323846

void velocity_rand(Mol &mol){
    mol.velocity.x = cos(2 * M_PI * static_cast<double>(rand()) / RAND_MAX);
    mol.velocity.y = cos(2 * M_PI * static_cast<double>(rand()) / RAND_MAX);
}

#endif