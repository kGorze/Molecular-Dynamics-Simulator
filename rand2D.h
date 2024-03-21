#ifndef rand2D_h
#define rand2D_h
#include <cmath>
#include <random>

struct Mol;

void myFunction(Mol *mol){
    mol->rv.x = std::cos(2*cmath::M_PI*rand()/RAND_MAX);
    mol->rv.y = std::cos(2*cmath::M_PI*rand()/RAND_MAX);
}

#endif