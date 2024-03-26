/*
This chapter provides the introductory appetizer and aims to leave the reader new to MD with a feeling for what the subject is all about. Later chapters will address the techniques in detail; here the goal is to demonstrate a working example with a minimum of fuss and so convince the beginner that MD is not only straightfor- ward but also that it works successfully. 

Of course, the technique for evaluating the forces discussed here is not particularly efficient from a computational point of view and the model is about the simplest there is
*/

#include "../MDS/hfiles/vector_operations.h"
#include "../MDS/hfiles/types.h"
#include "../MDS/hfiles/rand2D.h"
#include "../MDS/hfiles/statistics.h"
#include "../MDS/hfiles/physical.h"
#include "../MDS/hfiles/leapfrog.h"


#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

#define Sqr(x) ((x)*(x)) 
#define Cube(x) ((x)*(x)*(x)) //macros for squaring and cubing numbers
#define DO_MOL for(n = 0; n< nMol;n++) //This is a macro that is used to loop over all the molecules in the system. It is used in the ComputeForces and EvalProps functions.
#define NDIM 2

#define AllocMem(a,n,t) a = (t*)malloc((n)*sizeof(t)) //This is a macro that is used to allocate memory for arrays. It is used in the main program to allocate memory for the mol array.





#define VSCopy(v2,s1,v1){ (v2).x = (s1) * (v1).x; (v2).y = (s1) * (v1).y;}


/*
 From a practical point of view, writ- ing *mol in the above list of declarations is equivalent to mol[...] with a specific array size, except that in the former case the array size is established when the program is run rather than at compilation time
*/
VecR region, vSum,initUcell; 
/*
The vector region contains the edge lengths of the simulation region.
*/
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uCut, uSum, velMag, virSum, vvSum;
int moreCycels, stepAvg, stepCount, stepEquil, stepLimit;

int nMol = 400;
Mol* mol = (Mol*)malloc(sizeof(Mol)*nMol); // Change the variable mol to be a vector of Mols// Change the allocation of memory for mol to use the vector class


FILE* filePtr;

// void AllocArrays(){
//     AllocMem(mol,nMol,Mol);
// };


// This function calculates the interatomic forces between all pairs of atoms
// within a cutoff distance for a Molecular Dynamics simulation.

void ComputeForces(){
    VecR dr, fc; // Vector variables for distance between atoms and force between them
    real fcVal, rr, rrCut, rri, rri3; // Variables for force magnitude, distance squared, cutoff distance squared, and reciprocal distance terms
    int j1, j2, n; // Loop variables and atom index

    rrCut = Sqr(rCut); // Square of the cutoff distance for efficient comparison
    DO_MOL{VZero(&mol[n].accelaration);} // Zero out accelerations for all atoms
    uSum = 0; // Initialize potential energy sum
    virSum = 0; // Initialize virial sum

    // Nested loops over all pairs of atoms
    for(j1 = 0; j1 < nMol-1; j1++){
        for(j2 = j1+1; j2 < nMol; j2++){
            
            // Calculate distance vector between atom j1 and j2
            VSub(&dr, &mol[j1].coordinates, &mol[j2].coordinates);
            // Apply periodic boundary conditions to ensure the shortest distance
            VWrapAll(dr, region);
            // Calculate squared distance between atoms
            rr = VLenSq(&dr);
            // Check if the distance is within the cutoff radius
            if(rr < rrCut){
                // Calculate reciprocal distance terms
                rri = 1.0 / rr;
                rri3 = Cube(rri);
                // Calculate Lennard-Jones force magnitude
                fcVal = 48.0 * rri3 * (rri3 - 0.5) * rri;
                // Apply forces to atoms j1 and j2
                VVSAdd(&mol[j1].accelaration, fcVal, &dr);
                VVSAdd(&mol[j2].accelaration, -fcVal, &dr);
                // Accumulate potential energy and virial contribution
                uSum += 4.0 * rri3 * (rri3 - 1.0) - uCut;
                virSum += fcVal * rr;
            }
        }
    }
}



void LeapfrogStep(int part){
    int n;
    if(part == 1){
        DO_MOL{
            leapfrog_velocity(mol[n].velocity, 0.5*deltaT, mol[n].accelaration);
            leapfrog_coordinates(mol[n].coordinates, deltaT, mol[n].velocity);
        }
    }else{
        DO_MOL{
            leapfrog_velocity(mol[n].velocity, 0.5*deltaT, mol[n].accelaration);
        }
    }
};



void ApplyBoundaryCond(){
    //responsible for taking care of any periodic wraparound in the updated coordinates
    int n;
    DO_MOL{
        VWrapAll(mol[n].coordinates, region);
    }

};

void InitCoords(){
    //initializes the coordinates of the particles
    VecR c;
    int n, nx, ny;

    VecR gap;
    VDiv(&gap,&region,&initUcell);
    n= 0;
    for(nx = 0; nx<initUcell.x; nx++){
        for(ny = 0; ny<initUcell.y; ny++){
            
            VSet(&c, nx+0.5, ny+0.5);
            VMul(&c,&c,&gap);
            VVSAdd(&c, -0.5, &region);
            Mol temp = {c, {0,0}, {0,0}};
            mol[n] = temp;
            n++;
        }
    }
}


void InitVels(){
    //initializes the velocities of the particles
    int n;
    VZero (&vSum); //accumulate the total velocity(momentum)
    DO_MOL{
        velocity_rand(&mol[n]);
        VScale(&mol[n].velocity, velMag);
        VVAdd(&vSum, &mol[n].velocity);
    }

    DO_MOL{
        VVSAdd(&mol[n].velocity, (-1.0/nMol), &vSum);
    }
};

void InitAccels(){
    //initializes the accelerations of the particles
    int n;
    DO_MOL{
        VZero(&mol[n].accelaration);
    }
};


void EvalProps(){
    real vv;
    int n;

    VZero(&vSum);

    vvSum = 0;
    DO_MOL{
        VVAdd(&vSum, &mol[n].velocity);
        vv = VLenSq(&mol[n].velocity);
        vvSum += vv;
    }
    kinEnergy.val = (0.5*vvSum)/nMol;
    totEnergy.val = kinEnergy.val + (uSum/nMol);
    pressure.val = density*(vvSum - virSum)/(nMol*NDIM);
};


void AccumProps(int icode){
    if(icode == 0){
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
};

void SetParams(vector<KeyValue> *data){
    for (const auto &param : {"deltaT", "density", "initUcell", "stepAvg", "stepEquil", "stepLimit", "temperature"}) {
        // Find the entry in data that corresponds to the current parameter
        string key = param;
        string value;
        for (const auto &entry : (*data)) {
            if (entry.key == key) {
                size_t pos = entry.value.find_first_not_of(" ");
                value = entry.value.substr(pos);
                break;
            }
        }
        //cout<<key<<value;
        // Use the value to set the corresponding parameter
        if (key == "deltaT") {
            deltaT = stod(value);
        } else if (key == "density") {
            density = stod(value);
        } else if (key == "initUcell") {
            // Assuming initUcell is a vector of two integers
            int x, y;
            sscanf(value.c_str(), "%d %d", &x, &y);
            initUcell.x = x;
            initUcell.y = y;
        } else if (key == "stepAvg") {
            stepAvg = stoi(value);
        } else if (key == "stepEquil") {
            stepEquil = stoi(value);
        } else if (key == "stepLimit") {
            stepLimit = stoi(value);
        } else if (key == "temperature") {
            temperature = stod(value);
        }
    }

    rCut = pow(2.0, 1.0/6.0);

    VSCopy(region, 1/sqrt(density), initUcell);
    nMol = VProd(&initUcell); //The evaluation of nMol and region assumes just one atom per unit cell, and allowance is made for momentum conservation
    //(which removes NDIM degrees of freedom)
    velMag = sqrt(NDIM*(1 - (1/nMol)*temperature));
};

void PrintSummary (FILE *fp)
{

// fprintf(fp, "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
//         stepCount, timeNow, VCSum(&vSum) / nMol, PropEst(totEnergy),
//         PropEstSig(totEnergy), PropEst(kinEnergy), PropEstSig(kinEnergy),
//         PropEst(pressure), PropEstSig(pressure));
// fflush(fp);
// fflush(stdout);       


cout<<stepCount<<" "<<timeNow<<" "<<VCSum(&vSum) / nMol<<" "
    <<PropEst(totEnergy)<<" "<<PropEstSig(totEnergy)
    <<" "<<PropEst(kinEnergy)<<" "<<PropEstSig(kinEnergy)
    <<" "<<PropEst(pressure)<<" "<<PropEstSig(pressure)<<endl;
}; 

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
        PrintSummary(filePtr);
        AccumProps(0);
    }

};


void SetupJob(){
    stepCount = 0;
    InitCoords(); // initial coordinates
    InitVels(); // initial velocities
    InitAccels(); // initial accelerations
    AccumProps(0);
};



void GetNameList(const char* fd, vector<KeyValue>* data) {
    ifstream file(fd);
    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
        return;
    }

    string line;
    const string pattern = "initUcell";

    while (getline(file, line)) {
        KeyValue kv;

        if (line.substr(0, pattern.size()) == pattern) {
            line.erase(0, pattern.size());
            kv.key = pattern;
            size_t pos = line.find_first_not_of(" \t");
            if (pos != string::npos) {
                line.erase(0, pos);
                kv.value = line;
                (*data).push_back(kv);
            }
        } else {
            size_t pos = line.find_first_of(" \t");
            if (pos != string::npos) {
                kv.key = line.substr(0, pos);
                line.erase(0, pos);
                size_t pos = line.find_first_not_of(" ");
                line.erase(0, pos);
                kv.value = line;
                (*data).push_back(kv);
            }
        }
    }

    file.close();
};

void PrintNameList (vector<KeyValue> *data){
    for(const auto& kv : (*data)) {
        cout << kv.key << " " << kv.value << endl;
    }
};



int main(){

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
    
    FILE* filePtr = fopen("summary.txt", "w"); // Open file in write mode
    if (filePtr == nullptr) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return 1;
    }
    
    vector<KeyValue> data;
    GetNameList("data.in", &data);
    //PrintNameList(&data);
    SetParams(&data);
    InitCoords();
    SetupJob();

    //print all moll coordinates
    // for(int i = 0; i < mol.size(); i++){
    //     cout<<mol[i].coordinates.x<<" "<<mol[i].coordinates.y<<endl;
    // }
    
    int moreCycles = 1;
    while(moreCycles){
        SingleStep();
        if(stepCount>=stepLimit){
            moreCycles = 0; 
        }
    }
    free(mol);
    fclose(filePtr);
    return 0;
};