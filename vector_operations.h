#ifndef vector_operations_h
#define vector_operations_h

#define VAdd(v1,v2,v3){(v1).x = (v2).x + (v3).x; (v1).y = (v2).y + (v3).y;} //addition of two vectors
#define VSub(v1,v2,v3){(v1).x = (v2).x - (v3).x; (v1).y = (v2).y - (v3).y;} //subtraction of two vectors
#define VMul(v1,v2,v3){(v1).x = (v2).x*(v3).x;(v1).y = (v2).y*(v3).y;} //multiplication of two vectors
#define VDot(v1,v2) ((v1).x*(v2).x + (v1).y*(v2).y) //dot product of two vectors

#endif