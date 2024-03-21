#ifndef TYPES_H
#define TYPES_H
#include <string>

typedef struct{
    double x;
    double y;
} VecR;

typedef double real;


typedef struct {
    std::string key;
    std::string value;
} KeyValue;


typedef struct{
    VecR coordinates; //coordinates
    VecR velocity; //velocity
    VecR accelaration; //acceleration
} Mol;

typedef struct {
    int x,y;
} VecI;

typedef struct {
    double val, sum, sum2;
} Prop;

#endif