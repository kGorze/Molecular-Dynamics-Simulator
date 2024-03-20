/*
This chapter provides the introductory appetizer and aims to leave the reader new to MD with a feeling for what the subject is all about. Later chapters will address the techniques in detail; here the goal is to demonstrate a working example with a minimum of fuss and so convince the beginner that MD is not only straightfor- ward but also that it works successfully. 

Of course, the technique for evaluating the forces discussed here is not particularly efficient from a computational point of view and the model is about the simplest there is
*/

#include <iostream>
#include <math.h>
using namespace std;

#define Sqr(x) ((x)*(x))
#define Cube(x) ((x)*(x)*(x))
#define DO_MOL for(n = 0; n< nMol;n++)

#define VAdd(v1,v2,v3){(v1).x = (v2).x + (v3).x; (v1).y = (v2).y + (v3).y;}
#define VSub(v1,v2,v3){(v1).x = (v2).x - (v3).x; (v1).y = (v2).y - (v3).y;}
#define VDot(v1,v2) ((v1).x*(v2).x + (v1).y*(v2).y)
#define VSAdd(v1,v2,s3,v3){(v1).x = (v2).x + (s3)*(v3).x; (v1).y = (v2).y + (s3)*(v3).y;}
#define VSet(v1,s2){(v1).x = s2; (v1).y = s2;}
#define VSetAll(v,s)(VSET(v,s,s))
#define VZero(v)(VSetAll(v,0.0))
#define VVSAdd(v1,s2,v2)(VSAdd(v1,v1,s2,v2))
#define VLenSq(v)(VDot(v,v))

#define VWrap(v,t){
    if(v.t >= 0.5*region.t){
        v.t -= region.t;
    }else{
        if(v.t < -0.5*region.t){
            v.t += region.t;
        }
    }
}
#define VWrapAll(v){VWrap(v,x); VWrap(v,y);}


typedef double real;

typedef struct{
    real x, y, z;
} VecR;

typedef struct{
    VecR r; //coordinates
    VecR rv; //velocity
    VecR ra; //acceleration
} Mol;

// void ComputeForces(){
//     VecR dr;
//     real fcVal, rr, rrCut, rri, rri3;
//     int j1, j2, n;

//     rrCut = sqrt(rCut);
//     for(n = 0; n< nMol; n++){
//         mol[n].ra.x = 0.0;
//         mol[n].ra.y = 0.0;
//     }
//     uSum = 0;
//     for(j1 = 0; j1 < nMol-1; j1++){
//         for(j2 = j1+1; j2 < nMol; j2++){
//             dr.x = mol[j1].r.x - mol[j2].r.x;
//             dr.y = mol[j1].r.y - mol[j2].r.y;

//             //WRAP THE PARTICLES BACK INTO THE BOX
//             if(dr.x > 0.5*region.x){
//                 dr.x -= region.x;
//             }
//             if(dr.x < -0.5*region.x){
//                 dr.x += region.x;
//             }
//             if(dr.y > 0.5*region.y){
//                 dr.y -= region.y;
//             }
//             if(dr.y < -0.5*region.y){
//                 dr.y += region.y;
//             }
//             rr = dr.x*dr.x + dr.y*dr.y;
//             if(rr < rCut2){
//                 rri = 1.0/rr;
//                 rri3 = rri*rri*rri;
//                 fcVal = 48.0*rri3*(rri3-0.5)*rri;
//                 mol[j1].ra.x += fcVal*dr.x;
//                 mol[j1].ra.y += fcVal*dr.y;
//                 mol[j2].ra.x -= fcVal*dr.x;
//                 mol[j2].ra.y -= fcVal*dr.y;
//                 uSum += 4.0*rri3*(rri3-1.0) - uCut;
//             }
//         }
//     }
// }

void ComputeForces(){
    VecR dr;
    real fcVal, rr, rrCut, rri, rri3;
    int j1, j2, n;

    rrCut = Sqr(rCut);
    DO_MOL{VZero(mol[n].ra);}
    uSum = 0;
    virSum = 0;
    for(j1 = 0;j1< nMol-1;j1++){
        for(j2 = j1+1;j2<nMol;j2++){
            VSub(dr, mol[j1].r, mol[j2].r);
            VWrapAll(dr);
            rr = VLenSq(dr);
            if(rr < rrCut){
                rri = 1.0/rr;
                rri3 = Cube(rri);
                fcVal = 48.0*rri3*(rri3-0.5)*rri;
                VVSAdd(mol[j1].ra, fcVal, dr);
                VVSAdd(mol[j2].ra, -fcVal, dr);
                uSum += 4.0*rri3*(rri3-1.0) - uCut;
                virSum += fcVal*rr;
            }
        }
    }
};

void LeapfrogStep(int part){
    // handles the task of integrating the coordinates and velocities
    int n;
    real s;
    if(part == 1){
        s = 0.5*deltaT;
        DO_MOL{
            VSAdd(mol[n].rv, s, mol[n].ra);
            VSAdd(mol[n].r, deltaT, mol[n].rv);
        }
    }else{
        s = 0.5*deltaT;
        DO_MOL{
            VSAdd(mol[n].rv, s, mol[n].ra);
        }
    }
}

void ApplyBoundaryCond(){
    //responsible for taking care of any periodic wraparound in the updated coordinates
    int n;
    DO_MOL{
        VWrapAll(mol[n].r);
    }

}


void SingleStep(){
    stepCount++;
    timeNow = stepCount* deltaT;
    LeapfrogStep(1);
    ApplyBoundaryCond();
    
    ComputeForces();



    LeapfrogStep(2);
    EvalProps();
    AccumProps(1);
    if(stepCount % stepAvg == 0){
        AccumProps(2);
        PrintSummary();
        AccumProps(0);
    }

}

void StupJob(){
    AllocArrays();
    stepCount = 0;
    InitCords();
    InitVels();
    InitAccels();
    AccumProps(0);
}



int main(){
    GetNameList(argc, argv);
    PrintNameList(stdout);
    SetParams();
    SetupJob();
    int moreCycles = 1;
    while(moreCycles){
        SingleStep();
        if(stepCount>=stepLimit){
            moreCycles = 0;
        }
    }
}