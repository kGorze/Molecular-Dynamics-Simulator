#ifndef vector_operations_h
#define vector_operations_h

typedef struct{
    double x;
    double y;
} Vector;

void VAdd(Vector *v1, const Vector *v2, const Vector *v3);  //addition of two vectors
void VSub(Vector *v1, const Vector *v2, const Vector *v3);  //subtraction of two vectors
void VMul(Vector *v1, const Vector *v2, const Vector *v3);  //multiplication of two vectors

float VDot(Vector *v1, const Vector *v2);  //dot product of two vectors
void VSAdd(Vector *v1, const Vector *v2, double s3, const Vector *v3);  //scalar addition of two vectors

void VSet(Vector *v1, double s2);  //set the value of a vector
void VZero(Vector *v);  //set the value of all the components of a vector to zero

void VVSAdd(Vector *v1, double s2, const Vector *v2);  //add a scalar multiple of one vector to another
void VScale(Vector *v1, double s2);

void VScale(Vector *v1, double s2){
    v1->x = v1->x * s2;
    v1->y = v1->y * s2;
}

void VVSAdd(Vector *v1, double s2, const Vector *v2) {
    VSAdd(v1, v1, s2, v2);
}

void VZero(Vector *v){
    v->x = 0;
    v->y = 0;
}

void VSet(Vector *v1, double s2){
    v1->x = s2;
    v1->y = s2;
}

void VAdd(Vector *v1, const Vector *v2, const Vector *v3){
    v1->x = v2->x + v3->x;
    v1->y = v2->y + v3->y;
}

void VSub(Vector *v1, const Vector *v2, const Vector *v3){
    v1->x = v2->x - v3->x;
    v1->y = v2->y - v3->y;
}

void VMul(Vector *v1, const Vector *v2, const Vector *v3){
    v1->x = v2->x * v3->x;
    v1->y = v2->y * v3->y;
}

float VDot(Vector *v1, const Vector *v2){
    return v1->x * v2->x + v1->y * v2->y;
}

void VSAdd(Vector *v1, const Vector *v2, double s3, const Vector *v3){
    v1->x = v2->x + s3*v3->x; v1->y = v2->y + s3*v3->y;
}

#endif