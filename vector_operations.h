#pragma once
#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H
#include "types.h"


#define Sqr(x) ((x)*(x)) 
#define Cube(x) ((x)*(x)*(x)) //macros for squaring and cubing numbers
#define DO_MOL for(this->mnDOMOL = 0; this->mnDOMOL< this->mnMol;this->mnDOMOL++) //This is a macro that is used to loop over all the molecules in the system. It is used in the ComputeForces and EvalProps functions.
#define NDIM 2
#define Min(x,y) ((x)<(y)?(x):(y)) //This is a macro that is used to find the minimum of two numbers. It is used in the EvalVelDist function.
#define Max(x,y) ((x)>(y)?(x):(y)) //This is a macro that is used to find the maximum of two numbers. It is used in the EvalVelDist function.

//#define AllocMem(a,n,t) a = (t*)malloc((n)*sizeof(t))
//#define VSCopy(v2,s1,v1)((v2).x = (s1) * (v1).x (v2).y = (s1) * (v1).y)


void VAdd(VecR & v1, const VecR & v2, const VecR &v3);  //addition of two VecRs
void VSub(VecR & v1, const VecR & v2, const VecR &v3);  //subtraction of two VecRs
void VMul(VecR & v1, const VecR & v2, const VecR &v3);  //multiplication of two VecRs
double VLenSq(VecR & v1);  //length of a VecR squared

double VDot(VecR & v1, const VecR & v2);  //dot product of two VecRs
double VDotRef(VecR v1, VecR v2);  //dot product of two VecRs
void VSAdd(VecR & v1, const VecR & v2, double s3, const VecR & v3);  //scalar addition of two VecRs

void VSet(VecR & v1, double s2, double s3);  //set the value of a VecR
void VZero(VecR & v);  //set the value of all the components of a VecR to zero

void VVSAdd(VecR & v1, double s2, const VecR & v2);  //add a scalar multiple of one VecR to another
void VScale(VecR & v1, double s2);
void VSCopy(VecR & v1, double s2, const VecR & v2);
void VDiv(VecR & result, VecR & numerator, VecR & denominator);  //element-wise division of two VecRs
void VVAdd(VecR & v1, VecR & v2);  //element-wise addition of two VecRs
double VCSum(VecR & v1);  //sum of the components of a VecR
double VProd(VecR & v1); //product of the components of a VecR
double VLen(VecR & v1);  //length of a VecR

#endif