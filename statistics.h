#pragma once
#ifndef STATISTICS_H
#define STATISTICS_H

#include "vector_operations.h"

// Declaration to zero out the properties represented by v
void PropZero(Prop &v);

// Declaration to accumulate values into the properties represented by v
void PropAccum(Prop &v);

// Declaration to calculate the average and standard deviation of the properties represented by v
void PropAvg(Prop &v, int n);

// Declaration to compute an estimate of the property represented by v
double PropEst(const Prop &v);

// Declaration to compute an estimate of the standard deviation of the property represented by v
double PropEstSig(const Prop &v);

#endif