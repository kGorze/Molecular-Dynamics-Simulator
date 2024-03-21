#ifndef STATISTICS_H
#define STATISTICS_H

#include <cmath>
#include <algorithm>
#include "types.h"

// Define the structure for properties


// Function to zero out the properties represented by v
inline void PropZero(Prop &v) {
    v.sum = 0.0;
    v.sum2 = 0.0;
}

// Function to accumulate values into the properties represented by v
inline void PropAccum(Prop &v) {
    v.sum += v.val;
    v.sum2 += v.val * v.val;
}

// Function to calculate the average and standard deviation of the properties represented by v
inline void PropAvg(Prop &v, int n) {
    v.sum /= n;
    v.sum2 = sqrt(std::max(v.sum2 / n - v.sum * v.sum, 0.0));
}

// Function to compute an estimate of the property represented by v
inline double PropEst(const Prop &v) {
    return v.sum;
}

inline double PropEstSig(const Prop &v) {
    return v.sum2;
}

#endif // PROP_FUNCTIONS_H
