#ifndef vector_operations_h
#define vector_operations_h

typedef struct{
    double x;
    double y;
} VecR;

void VAdd(VecR *v1, const VecR *v2, const VecR *v3);  //addition of two VecRs
void VSub(VecR *v1, const VecR *v2, const VecR *v3);  //subtraction of two VecRs
void VMul(VecR *v1, const VecR *v2, const VecR *v3);  //multiplication of two VecRs

float VDot(VecR *v1, const VecR *v2);  //dot product of two VecRs
void VSAdd(VecR *v1, const VecR *v2, double s3, const VecR *v3);  //scalar addition of two VecRs

void VSet(VecR *v1, double s2, double s3);  //set the value of a VecR
void VZero(VecR *v);  //set the value of all the components of a VecR to zero

void VVSAdd(VecR *v1, double s2, const VecR *v2);  //add a scalar multiple of one VecR to another
void VScale(VecR *v1, double s2);
void VSCopy(VecR *v1, double s2, const VecR *v2);
void VDiv(VecR *result, VecR *numerator,  VecR *denominator);  //element-wise division of two VecRs

int VProd(VecR *v1);

int VProd(VecR *v1){
    return v1->x * v1->y;
}

void VDiv(VecR *result, VecR *numerator, VecR *denominator){
    result->x = numerator->x / denominator->x;
    result->y = numerator->y / denominator->y;
}

void VSCopy(VecR *v1, double s2, const VecR *v2){
    v1->x = s2 * v2->x;
    v1->y = s2 * v2->y;
}


void VScale(VecR *v1, double s2){
    v1->x = v1->x * s2;
    v1->y = v1->y * s2;
}

void VVSAdd(VecR *v1, double s2, const VecR *v2) {
    VSAdd(v1, v1, s2, v2);
}

void VZero(VecR *v){
    v->x = 0;
    v->y = 0;
}

void VSet(VecR *v1, double s2, double s3){
    v1->x = s2;
    v1->y = s3;
}

void VAdd(VecR *v1, const VecR *v2, const VecR *v3){
    v1->x = v2->x + v3->x;
    v1->y = v2->y + v3->y;
}

void VSub(VecR *v1, const VecR *v2, const VecR *v3){
    v1->x = v2->x - v3->x;
    v1->y = v2->y - v3->y;
}

void VMul(VecR *v1, const VecR *v2, const VecR *v3){
    v1->x = v2->x * v3->x;
    v1->y = v2->y * v3->y;
}

float VDot(VecR *v1, const VecR *v2){
    return v1->x * v2->x + v1->y * v2->y;
}

void VSAdd(VecR *v1, const VecR *v2, double s3, const VecR *v3){
    v1->x = v2->x + s3*v3->x; v1->y = v2->y + s3*v3->y;
}

#endif