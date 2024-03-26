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

real rangeVel;
int countVel, limitVel, stepVel;

int nMol = 400;
int sizeHistVel = 50;

vector<Mol> mol(nMol);
vector<real> histVel(sizeHistVel);


FILE* filePtr;
FILE* histoPtr;




// void AllocArrays(){
//     AllocMem(mol,nMol,Mol);
// };


// This function calculates the interatomic forces between all pairs of atoms
// within a cutoff distance for a Molecular Dynamics simulation.

void EvalVelDist();
void ComputeForces();
void LeapfrogStep(int part);
void ApplyBoundaryCond();
void InitCoords();
void InitVels();
void InitAccels();
void EvalProps();
void AccumProps(int icode);
void SetParams(vector<KeyValue>* data);
void PrintSummary(FILE *fp);
void SetupJob();
void SingleStep(FILE *fp);
void GetNameList(const char *fileName, vector<KeyValue> *data);
void PrintNameList(vector<KeyValue> *data);
void EvalVelDist();
void PrintVelDist(FILE *fp);






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
        std::cerr << "Error opening file for writing summary" << std::endl;
        return 1;
    }

    FILE* histoPtr = fopen("histo.txt", "w"); // Open file in write mode
    if (histoPtr == nullptr) {
        std::cerr << "Error opening file for writing! histo" << std::endl;
        return 1;
    }


    vector<KeyValue> data;
    GetNameList("data.in", &data);
    
    //PrintNameList(&data);
    SetParams(&data);

    
    cout<<"test"<<endl;


    InitCoords();
    SetupJob();


    //print coordinates of the molecules
    // for(int i = 0; i < nMol; i++){
    //     cout<<mol[i].coordinates.x<<" "<<mol[i].coordinates.y<<endl;
    // }
    
    int moreCycles = 1;
    while(moreCycles){
        SingleStep(filePtr);
        //temperature +=0,5;
        if(stepCount>=stepLimit){
            moreCycles = 0; 
        }
    }
    PrintVelDist(histoPtr);

    fclose(filePtr);
    fclose(histoPtr);
    return 0;
};

void ComputeForces ()
{
    VecR dr; // Vector representing the distance between two atoms
    real fcVal, rr, rrCut, rri, rri3; // Variables for force magnitude, distance squared, cutoff distance squared, and reciprocal distance terms
    int j1, j2, n; // Loop variables and atom index

    rrCut = Sqr (rCut); // Calculate the square of the cutoff distance for efficient comparison

    // Initialize the acceleration components of each atom to zero
    for (n = 0; n < nMol; n++) {
        mol[n].accelaration.x = 0.;
        mol[n].accelaration.y = 0.;
    }

    uSum = 0.; // Initialize the potential energy sum

    // Nested loops over all pairs of atoms
    for (j1 = 0; j1 < nMol - 1; j1++) {
        for (j2 = j1 + 1; j2 < nMol; j2++) {
            // Calculate the distance vector between atom j1 and j2
            dr.x = mol[j1].coordinates.x - mol[j2].coordinates.x;
            dr.y = mol[j1].coordinates.y - mol[j2].coordinates.y;

            // Apply periodic boundary conditions to ensure the shortest distance
            if (dr.x >= 0.5 * region.x) dr.x -= region.x;
            else if (dr.x < -0.5 * region.x) dr.x += region.x;
            if (dr.y >= 0.5 * region.y) dr.y -= region.y;
            else if (dr.y < -0.5 * region.y) dr.y += region.y;

            // Calculate squared distance between atoms
            rr = dr.x * dr.x + dr.y * dr.y;

            // Check if the distance is within the cutoff radius
            if (rr < rrCut) {
                // Calculate reciprocal distance terms
                rri = 1. / rr;
                rri3 = rri * rri * rri;
                // Calculate Lennard-Jones force magnitude
                fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                // Apply forces to atoms j1 and j2
                mol[j1].accelaration.x += fcVal * dr.x;
                mol[j1].accelaration.y += fcVal * dr.y;
                mol[j2].accelaration.x -= fcVal * dr.x;
                mol[j2].accelaration.y -= fcVal * dr.y;
                // Accumulate potential energy
                uSum += 4. * rri3 * (rri3 - 1.) + 1.;
            }
        }
    }
}

//GOOD
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

//GOOD
void ApplyBoundaryCond(){
    //responsible for taking care of any periodic wraparound in the updated coordinates
    int n;
    DO_MOL{
        VWrapAll(mol[n].coordinates, region);
    }

};

//GOOD
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
            //Mol temp = {c, {0,0}, {0,0}};
            mol[n].coordinates = c;
            ++n;
        }
    }
}

//GOOD
void InitVels(){
    //initializes the velocities of the particles
    int n;
    VZero (&vSum); //accumulate the total velocity(momentum)

    DO_MOL{
        velocity_rand(&mol[n]);
        //cout<<"Velocity:";
        //cout<<mol[n].velocity.x<<" "<<mol[n].velocity.y<<endl;
        VScale(&mol[n].velocity, velMag);
        VVAdd(&vSum, &mol[n].velocity);
        //cout<<"Sum:";
        //cout<<vSum.x<<" "<<vSum.y<<endl;
    }

    DO_MOL{
        VVSAdd(&mol[n].velocity, (-1.0/nMol), &vSum);
    }

    // DO_MOL{
    //     cout<<"Velocity of atom nr:"<<n<<endl;
    //     cout<<mol[n].velocity.x<<" "<<mol[n].velocity.y<<endl;
    // }
    //cout<<"Sum:";
    //cout<<vSum.x<<" "<<vSum.y<<endl;
};

//GOOD
void InitAccels(){
    //initializes the accelerations of the particles
    int n;
    DO_MOL{
        VZero(&mol[n].accelaration);
    }
    // DO_MOL{
    //     cout<<"Acceleration of atom nr:"<<n<<endl;
    //     cout<<mol[n].accelaration.x<<" "<<mol[n].accelaration.y<<endl;
    // }
};

//GOOD
void EvalProps(){
     real vv; // Variable to store the squared velocity magnitude of each atom
    int n; // Loop variable representing the atom index

    VZero(&vSum); // Initialize the sum of velocities to zero

    vvSum = 0; // Initialize the sum of squared velocities to zero

    // Loop over all atoms in the system
    DO_MOL{
        // Accumulate the velocity vector of each atom
        VVAdd(&vSum, &mol[n].velocity);
        
        // Calculate the squared velocity magnitude of each atom
        vv = VLenSq(&mol[n].velocity);
        
        // Accumulate the squared velocity magnitude
        vvSum += vv;
    }

    // Calculate the kinetic energy per particle
    kinEnergy.val = (0.5 * vvSum) / nMol;

    // Calculate the total energy per particle
    totEnergy.val = kinEnergy.val + (uSum / nMol);

    // Calculate the pressure
    pressure.val = density * (vvSum - virSum) / (nMol * NDIM);
};

//GOOD
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

//GOOD
void SetParams(vector<KeyValue> *data){
    for (const auto &param : {"deltaT", "density", "initUcell", "stepAvg", "stepEquil", "stepLimit", "temperature", "limitVel", "rangeVel", "sizeHistVel", "stepVel"}) {
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
        } else if (key == "limitVel") {
            limitVel = stoi(value);
        } else if (key == "rangeVel") {
            rangeVel = stod(value);
        } else if (key == "sizeHistVel") {
            sizeHistVel = stoi(value);
        } else if (key == "stepVel") {
            stepVel = stoi(value);
        }
    }

    rCut = pow(2.0, 1.0/6.0);

    VSCopy(region, 1/sqrt(density), initUcell);
    nMol = VProd(&initUcell); //The evaluation of nMol and region assumes just one atom per unit cell, and allowance is made for momentum conservation
    //(which removes NDIM degrees of freedom)
    velMag = sqrt(NDIM*(1 - (1/nMol)*temperature));
}


void PrintSummary (FILE *fp)
{
    if(fp == nullptr){
        cerr<<"Error opening file"<<endl;
        return;
    }
    real vSumValue = VCSum(&vSum) / nMol;
    real totEnergyEst = PropEst(totEnergy);
    real totEnergySig = PropEstSig(totEnergy);
    real kinEnergyEst = PropEst(kinEnergy);
    real kinEnergySig = PropEstSig(kinEnergy);
    real pressureEst = PropEst(pressure);
    real pressureSig = PropEstSig(pressure);

    // Print the values using fprintf
    fprintf(fp, "%5d %8.4f %7.20f %7.20f %7.20f %7.20f %7.20f %7.20f %7.20f\n",
            stepCount, timeNow, vSumValue, totEnergyEst,
            totEnergySig, kinEnergyEst, kinEnergySig,
            pressureEst, pressureSig); 

cout<<stepCount<<" "<<timeNow<<" "<<VCSum(&vSum) / nMol<<" "
    <<PropEst(totEnergy)<<" "<<PropEstSig(totEnergy)
    <<" "<<PropEst(kinEnergy)<<" "<<PropEstSig(kinEnergy)
    <<" "<<PropEst(pressure)<<" "<<PropEstSig(pressure)<<endl;
 }; 

void SingleStep(FILE *fp){
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
        if(stepCount>=stepEquil && (stepCount - stepEquil) % stepVel == 0){
            EvalVelDist();
        }
        PrintSummary(fp);
        AccumProps(0);
    }

};


void SetupJob(){
    stepCount = 0;
    countVel = 0;
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

void EvalVelDist(){
    real deltaV, histSum;
    int j,n;

    if(countVel == 0){
        for(j = 0;j<sizeHistVel;j++){
            histVel[j] = 0;
        }
    }
    deltaV = rangeVel/sizeHistVel;
    DO_MOL{
        j = (VLen(&mol[n].velocity))/deltaV;
        ++histVel[Min(j,sizeHistVel-1)];
    }
    ++countVel;
    if(countVel == limitVel){
        histSum = 0;
        for(j = 0; j<sizeHistVel; j++){
            histSum += histVel[j];
        }
        for(j = 0; j<sizeHistVel; j++){
            histVel[j] /= histSum;
        }
        countVel = 0;
    }
    cout<<".";
}

void PrintVelDist(FILE *fp){
    real vBin;
    int n;

    printf ("vdist (%.3f)\n", timeNow);
    for (n = 0; n < sizeHistVel; n ++) {
        vBin = (n + 0.5) * rangeVel / sizeHistVel;
        fprintf (fp, "%8.3f %8.3f\n", vBin, histVel[n]);
    }
}