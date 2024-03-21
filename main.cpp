/*
This chapter provides the introductory appetizer and aims to leave the reader new to MD with a feeling for what the subject is all about. Later chapters will address the techniques in detail; here the goal is to demonstrate a working example with a minimum of fuss and so convince the beginner that MD is not only straightfor- ward but also that it works successfully. 

Of course, the technique for evaluating the forces discussed here is not particularly efficient from a computational point of view and the model is about the simplest there is
*/
#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;

#define Sqr(x) ((x)*(x)) 
#define Cube(x) ((x)*(x)*(x)) //macros for squaring and cubing numbers

#define DO_MOL for(n = 0; n< nMol;n++) //This is a macro that is used to loop over all the molecules in the system. It is used in the ComputeForces and EvalProps functions.

#define VAdd(v1,v2,v3){(v1).x = (v2).x + (v3).x; (v1).y = (v2).y + (v3).y;} //addition of two vectors
#define VSub(v1,v2,v3){(v1).x = (v2).x - (v3).x; (v1).y = (v2).y - (v3).y;} //subtraction of two vectors
#define VMul(v1,v2,v3){(v1).x = (v2).x*(v3).x;(v1).y = (v2).y*(v3).y;} //multiplication of two vectors
#define VDot(v1,v2) ((v1).x*(v2).x + (v1).y*(v2).y) //dot product of two vectors



#define VSAdd(v1,v2,s3,v3){(v1).x = (v2).x + (s3)*(v3).x; (v1).y = (v2).y + (s3)*(v3).y;} //scalar addition of two vectors


#define VSet(v1,s2){(v1).x = s2; (v1).y = s2;} //set the value of a vector


#define VSetAll(v,s)(VSet(v,s)) //set the value of all the components of a vector to the same value

#define VZero(v)(VSetAll(v,0)) //set the value of all the components of a vector to zero

#define VVSAdd(v1,s2,v2)(VSAdd(v1,v1,s2,v2)) //add a scalar multiple of one vector to another
#define VLenSq(v)(VDot(v,v)) //return the square of the length of a vector


//DO OGARNIECIA
#define VWrap(v,t){ if(v.t >= 0.5*region.t){ v.t -= region.t; else{ if(v.t < -0.5*region.t){ v.t += region.t;}}} 

#define VWrapAll(v){VWrap(v,x); VWrap(v,y);}

#define VScale(v,s){(v).x *= s;(v).y *= s;}
#define VVAdd(v1,v2) { VAdd(v1,v1,v2);}

#define NDIM 2

#define AllocMem(a,n,t) a = (t*)malloc((n)*sizeof(t))

#define VSCopy(v2,s1,v1){(v2).x = (s1)*(v1).x; (v2).y = (s1)*(v1).y;}
#define VProd(v) {(v).x * (v).y}

#define PropZero(v){v.sum = 0.0; v.sum2 = 0.0;}
#define PropAccum(v){v.sum += v.val; v.sum2 += Sqr(v.val);}
#define PropAvg(v,n){v.sum = /n; v.sum2 = sqrt(Max(v.sum2/n-Sqr(v.sum),0.0));}
#define PropEst(v) (v.sum/n)


typedef double real;

typedef struct{
    real x, y, z;
} VecR;

typedef struct{
    VecR r; //coordinates
    VecR rv; //velocity
    VecR ra; //acceleration
} Mol;

typedef struct {
    int x,y;
} VecI;

#define VSCopy(v2,s1,v1){ (v2).x = (s1) * (v1).x; (v2).y = (s1) * (v1).y;}
#define VProd{(v).x * (v).y}

typedef struct{
    real val, sum, sum2;
} Prop;

Mol *mol; // The variable mol is actually a pointer to a one-dimensional array that is allocated dynamically at the start of the run and sized according to the value of nMol. 
/*
 From a practical point of view, writ- ing *mol in the above list of declarations is equivalent to mol[...] with a specific array size, except that in the former case the array size is established when the program is run rather than at compilation time
*/
VecR region, vSum; 
/*
The vector region contains the edge lengths of the simulation region.
*/
VecI initUcell;
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uCut, uSum, virSum, velMag, virSum, vvSum;
int moreCycels, nMol, stepAvg, stepCount, stepEquil, stepLimit;


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

void InitCoords(){
    //initializes the coordinates of the particles
    VecR c, gap;
    int n, nx, ny;

    VDiv(gap,region,initUcell);
    n= 0;
    for(nx = 0; nx<initUcell.x; nx++){
        for(ny = 0; ny<initUcell.y; ny++){
    VSet(c, nx+0.5, ny+0.5);
    VMul(c,c,gap);
    VVSAdd(c, -0.5, region);
    mol[n].r = c;
    n++;
        }
    }
}


void InitVels(){
    //initializes the velocities of the particles
    int n;

    
    VZero (vSum); //accumulate the total velocity(momentum)
    DO_MOL{
        VRand(&mol[n].rv);
        VScale(mol[n].rv, velMag);
        VVAdd(vSum, mol[n].rv);

}
DO_MOL VVSAdd(mol[n].rv, -1.0/nMol, vSum);
}

void InitAccels(){
    //initializes the accelerations of the particles
    int n;
    DO_MOL VZero(mol[n].ra);
}

void AllocArrays(){
    AllocMem(mol,nMol,Mol);
}

void EvalProps(){
    real vv;
    int n;

    VZero(vSum);

    vvSum = 0;
    DO_MOL{
        VVAdd(vSum, mol[n].rv);
        vv = VLenSq(mol[n].rv);
        vvSum += vv;
    }
    kinEnergy.val = 0.5*vvSum/nMol;
    totEnergy.val = kinEnergy.val + uSum/nMol;
    pressure.val = density*(vvSum - virSum)/(nMol*NDIM);
}


void AccumProps(int icode){
    if(icode = 0){
        PropZero(totEnergy);
        PropZero(kinEnergy);
        PropZero(pressure);
    }
    else if (icode == 1){
        PropAccum(totEnergy);
        PropAccum(kinEnergy);
        PropAccum(pressure);
    }else if(icode == 2){
        PropAvg(totEnergy, stepAvg);
        PropAvg(kinEnergy, stepAvg);
        PropAvg(pressure, stepAvg);
    }
}

void SetParams(){
    rCut = pow(2.0,1.0/6.0); 
    VSCopy(region, 1/sqrt(density), initUcell);
    nMol = VProd(initUcell); //The evaluation of nMol and region assumes just one atom per unit cell, and allowance is made for momentum conservation
    //(which removes NDIM degrees of freedom)
    velMag = sqrt(NDIM*(1 - 1/nMol)*temperature);
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


void SetupJob(){
    AllocArrays();
    stepCount = 0;
    InitCords(); // initial coordinates
    InitVels(); // initial velocities
    InitAccels(); // initial accelerations
    AccumProps(0);
}

void PrintSummary (FILE *fp)
{
fprintf (fp,
"%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
stepCount, timeNow, VCSum (vSum) / nMol, PropEst (totEnergy),
PropEst (kinEnergy), PropEst (pressure));
}


int main(){

    real deltaT, density, temperature;
    int initUcell, stepAvg, stepEquil, stepLimit;

    /*
    deltaT 0.005
    density 0.8
    initUcell 20 20
    stepAvg 100
    stepEquil 0
    stepLimit 10000
    temperature 1    
    */

   //The initial configuration is a 20 × 20 square lattice so that there are a total of 400 atoms.
   //The timestep value deltaT is determined by the requirement that energy be conserved by the leapfrog method 
   //Temperature 1

//set parameters from input to the program
    GetNameList(argc, argv);
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