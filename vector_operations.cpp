#include "vector_operations.h"

double VCSum(VecR& v1) {
    return v1.x + v1.y;

}

double VLenSq(VecR& v1) {
    return VDotRef(v1, v1);
}

double VLen(VecR& v1) {
    return Sqr(VDotRef(v1, v1));
}


void VVAdd(VecR& v1, VecR& v2) {
    VAdd(v1, v1, v2);
}

double VProd(VecR& v1) {
    return v1.x * v1.y;
}

void VDiv(VecR& result, VecR& numerator, VecR& denominator) {
    result.x = numerator.x / denominator.x;
    result.y = numerator.y / denominator.y;
}

void VSCopy(VecR& v1, double s2, const VecR& v2) {
    v1.x = s2 * v2.x;
    v1.y = s2 * v2.y;
}


void VScale(VecR &v1, double s2) {
    v1.x *= s2;
    v1.y *= s2;
}

void VVSAdd(VecR &v1, double s2, const VecR& v2) {
    VSAdd(v1, v1, s2, v2);
}

void VZero(VecR& v) {
    v.x = 0;
    v.y = 0;
}

void VSet(VecR& v1, double s2, double s3) {
    v1.x = s2;
    v1.y = s3;
}

void VAdd(VecR& v1, const VecR& v2, const VecR& v3) {
    v1.x = v2.x + v3.x;
    v1.y = v2.y + v3.y;
}

void VSub(VecR& v1, const VecR& v2, const VecR& v3) {
    v1.x = v2.x - v3.x;
    v1.y = v2.y - v3.y;
}

void VMul(VecR &v1, const VecR& v2, const VecR& v3) {
    v1.x = v2.x * v3.x;
    v1.y = v2.y * v3.y;
}

double VDot(VecR& v1, VecR& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

double VDotRef(VecR v1, VecR v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

void VSAdd(VecR& v1, const VecR& v2, double s3, const VecR& v3) {
    v1.x = v2.x + s3 * v3.x;
    v1.y = v2.y + s3 * v3.y;
}
