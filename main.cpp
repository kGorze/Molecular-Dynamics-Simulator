/*
This chapter provides the introductory appetizer and aims to leave the reader new to MD with a feeling for what the subject is all about. Later chapters will address the techniques in detail; here the goal is to demonstrate a working example with a minimum of fuss and so convince the beginner that MD is not only straightfor- ward but also that it works successfully. 

Of course, the technique for evaluating the forces discussed here is not particularly efficient from a computational point of view and the model is about the simplest there is
*/

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef IOSTREAM_H
#define IOSTREAM_H
#include <iostream>
#endif

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef FSTREAM_H
#define FSTREAM_H
#include <fstream>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef UNORDERED_MAP_H
#define UNORDERED_MAP_H
#include <unordered_map>
#endif

#ifndef SSTREAM_H
#define SSTREAM_H
#include <sstream>
#endif

#ifndef CTIME_H
#define CTIME_H
#include <ctime>
#endif

#ifndef CSTDIO_H
#define CSTDIO_H
#include <cstdio>
#endif

#ifndef THREAD_H
#define THREAD_H
#include <thread>
#endif

#ifndef CHRONO_H
#define CHRONO_H
#include <chrono>
#endif

#include "inputhandler.h"
#include "vector_operations.h"
#include "types.h"
#include "rand2D.h"
#include "statistics.h"
#include "physical.h"
#include "leapfrog.h"
#include "render.h"

#define AllocMem(a,n,t) a = (t*)malloc((n)*sizeof(t)) //This is a macro that is used to allocate memory for arrays. It is used in the main program to allocate memory for the mol array.

#define VSCopy(v2,s1,v1){ (v2).x = (s1) * (v1).x; (v2).y = (s1) * (v1).y;}


/*
 From a practical point of view, writ- ing *mol in the above list of declarations is equivalent to mol[...] with a specific array size, except that in the former case the array size is established when the program is run rather than at compilation time
*/
VecR region, vSum,initUcell; 
/*
The std::vector region contains the edge lengths of the simulation region.
*/
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uCut, uSum, velMag, virSum, vvSum;
int moreCycels, stepAvg, stepCount, stepEquil, stepLimit;

real rangeVel;
int countVel, limitVel, stepVel;

int nMol = 400;
int sizeHistVel = 50;

std::vector<Mol> mol(nMol);
std::vector<real> histVel(sizeHistVel);

std::vector<std::vector<std::vector<double>>> dataCoords(500);

FILE* filePtr;
FILE* histoPtr;

real hFunction;


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
void SetParams(std::vector<KeyValue>* data);
void PrintSummary(FILE *fp);
void PrintCoordinates(const char* fileName);
std::vector<std::vector<std::vector<double>>> readFile(const std::string& filename);
void SetupJob();
double SingleStep(FILE *fp);
void GetNameList(const char *fileName, std::vector<KeyValue> *data);
void PrintNameList(std::vector<KeyValue> *data);
void EvalVelDist();
void PrintVelDist(FILE *fp);
int Simulation(unsigned int option);
void Showcase();



int main() {
    //Simulation(1); //if parameter is 1, the program will print the coordinates to a file
    //Showcase();
    InputHandler temp("data.in");
    return 0;
    };

void Showcase() {
    dataCoords = readFile("coordinates.txt");
    GLFWwindow* window = initializeScreen(800, 800);
    renderAtoms(window, dataCoords);
    terminateScreen();
}

int Simulation(unsigned int option){
    std::cout<<"-----------------------------------------------------------------------\n";
    std::cout <<"Initializing simulation\n";
    std::cout << "-----------------------------------------------------------------------\n\n";


    std::vector<KeyValue> data;

     //set parameters from input to the program
    FILE* filePtr = nullptr;
    FILE* histoPtr = nullptr;
    errno_t err;


    err = fopen_s(&filePtr, "summary.txt", "w"); // Open file in write mode
    if (err != 0) {
        std::cerr << "Error opening file for writing summary" << "\n";
        return 1;
    }

    err = fopen_s(&histoPtr, "histo.txt", "w"); // Open file in write mode
    if (err != 0) {
        std::cerr << "Error opening file for writing! histo" << "\n";
        return 1;
    }

    
    GetNameList("data.in", &data);
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "Parameters\n";
    std::cout << "-----------------------------------------------------------------------\n";
    PrintNameList(&data);
    std::cout << "-----------------------------------------------------------------------\n\n";
    
    std::this_thread::sleep_for(std::chrono::seconds(4));
    

    SetParams(&data);

    InitCoords();
    SetupJob();

    system("cls");

    int moreCycles = 1;
    float progress = 0.0;
    float stepFraction = 0.0;
    int barwidth = 50;

    while (moreCycles && progress < 1.0) {
        
        std::cout<< "Progress: [";
        int pos = barwidth * progress;
        for (int i = 0; i < barwidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
        std::cout<< "] " << int(progress * 100.0) << " %\r\n";
        std::flush(std::cout);
        SingleStep(filePtr);

        stepFraction = float(stepCount) / float(stepLimit);
        progress = stepFraction;

        if (stepCount >= stepLimit) {
            moreCycles = 0;
        }
    }
    //PrintVelDist(histoPtr);

    if (histoPtr != nullptr) {
        fclose(histoPtr);
    }
    if (filePtr != nullptr) {
        fclose(filePtr);
    }
    
    if (option == 1) {
		PrintCoordinates("coordinates.txt");
        return 0;
    }
    else {
        return 0;
    }
}

//make a function PrintCoordinates that will take name of the file - not a pointer to a file - to save the coordinates and open the file in write mode
void PrintCoordinates(const char* fileName) {
	FILE* filePtr = nullptr;
	errno_t err;

	err = fopen_s(&filePtr, fileName, "w"); // Open file in write mode
    if (err != 0) {
		std::cerr << "Error opening file for writing" << "\n";
		return;
	}

    //print the coordinates from dataCoords
    for (int i = 0; i < dataCoords.size(); i++) {
        for (int j = 0; j < dataCoords[i].size(); j++) {
			fprintf(filePtr, "%5d %8.4f %8.4f\n", i, dataCoords[i][j][0], dataCoords[i][j][1]);
		}
	}

    if (filePtr != nullptr) {
		fclose(filePtr);
	}

}

std::vector<std::vector<std::vector<double>>> readFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return {};
    }

    std::vector<std::vector<std::vector<double>>> dataCoords;

    std::string line;
    while (getline(file, line)) {
        std::istringstream iss(line);
        int index;
        double x, y;
        if (!(iss >> index >> x >> y)) {
            std::cerr << "Error parsing line: " << line << "\n";
            continue;
        }

        if (index >= dataCoords.size()) {
            dataCoords.resize(index + 1);
        }

        dataCoords[index].push_back({ x, y });
    }

    return dataCoords;
}




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
            // Calculate the distance std::vector between atom j1 and j2
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
            leapfrog_velocity(&mol[n].velocity, 0.5*deltaT, &mol[n].accelaration);
            leapfrog_coordinates(&mol[n].coordinates, deltaT, &mol[n].velocity);
        }
    }else{
        DO_MOL{
            leapfrog_velocity(&mol[n].velocity, 0.5*deltaT, &mol[n].accelaration);
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
        //std::cout<<"Velocity:";
        //std::cout<<mol[n].velocity.x<<" "<<mol[n].velocity.y<<"\n";
        VScale(&mol[n].velocity, velMag);
        VVAdd(&vSum, &mol[n].velocity);
        //std::cout<<"Sum:";
        //std::cout<<vSum.x<<" "<<vSum.y<<"\n";
    }

    DO_MOL{
        VVSAdd(&mol[n].velocity, (-1.0/nMol), &vSum);
    }

    // DO_MOL{
    //     std::cout<<"Velocity of atom nr:"<<n<<"\n";
    //     std::cout<<mol[n].velocity.x<<" "<<mol[n].velocity.y<<"\n";
    // }
    //std::cout<<"Sum:";
    //std::cout<<vSum.x<<" "<<vSum.y<<"\n";
};

//GOOD
void InitAccels(){
    //initializes the accelerations of the particles
    int n;
    DO_MOL{
        VZero(&mol[n].accelaration);
    }
    // DO_MOL{
    //     std::cout<<"Acceleration of atom nr:"<<n<<"\n";
    //     std::cout<<mol[n].accelaration.x<<" "<<mol[n].accelaration.y<<"\n";
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
        // Accumulate the velocity std::vector of each atom
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
void SetParams(std::vector<KeyValue> *data){
    std::unordered_map<std::string, std::string> paramMap;

    // Create a map of parameter names to values
    for (const auto &entry : *data) {
        paramMap[entry.key] = entry.value;
    }

    // Use the map to set the parameters
    if (paramMap.find("deltaT") != paramMap.end()) {
        deltaT = stod(paramMap["deltaT"]);
    }
    if (paramMap.find("density") != paramMap.end()) {
        density = stod(paramMap["density"]);
    }
    if (paramMap.find("initUcell") != paramMap.end()) {
        // Assuming initUcell is a std::vector of two integers
        int x, y;
        sscanf_s(paramMap["initUcell"].c_str(), "%d %d", &x, &y);
        initUcell.x = x;
        initUcell.y = y;
    }
    if (paramMap.find("stepAvg") != paramMap.end()) {
        stepAvg = stoi(paramMap["stepAvg"]);
    }
    if (paramMap.find("stepEquil") != paramMap.end()) {
        stepEquil = stoi(paramMap["stepEquil"]);
    }
    if (paramMap.find("stepLimit") != paramMap.end()) {
        stepLimit = stoi(paramMap["stepLimit"]);
    }
    if (paramMap.find("temperature") != paramMap.end()) {
        temperature = stod(paramMap["temperature"]);
    }
    if (paramMap.find("limitVel") != paramMap.end()) {
        limitVel = stoi(paramMap["limitVel"]);
    }
    if (paramMap.find("rangeVel") != paramMap.end()) {
        rangeVel = stod(paramMap["rangeVel"]);
    }
    if (paramMap.find("sizeHistVel") != paramMap.end()) {
        sizeHistVel = stoi(paramMap["sizeHistVel"]);
    }
    if (paramMap.find("stepVel") != paramMap.end()) {
        stepVel = stoi(paramMap["stepVel"]);
    }

    // Assuming rCut, region, nMol, and velMag are global variables
    rCut = pow(2.0, 1.0/6.0);
    VSCopy(region, 1/sqrt(density), initUcell);
    nMol = VProd(&initUcell);
    velMag = sqrt(NDIM*(1 - (1/nMol)*temperature));
}



void PrintSummary (FILE *fp)
{
    if(fp == nullptr){
        std::cerr<<"Error opening file"<<"\n";
        return;
    }
    real vSumValue = VCSum(&vSum) / nMol;
    real totEnergyEst = PropEst(totEnergy);
    real totEnergySig = PropEstSig(totEnergy);
    real kinEnergyEst = PropEst(kinEnergy);
    real kinEnergySig = PropEstSig(kinEnergy);
    real pressureEst = PropEst(pressure);
    real pressureSig = PropEstSig(pressure);

    fprintf(fp, " Step   Time    Sum(v)  Etot            EtotSig            Ekin            EkinSig            Pressure        PressureSig\n");
    fprintf(fp, "%5d %8.4f %7.20f %7.20f %7.20f %7.20f %7.20f %7.20f %7.20f\n",
            stepCount, timeNow, vSumValue, totEnergyEst,
            totEnergySig, kinEnergyEst, kinEnergySig,
            pressureEst, pressureSig); 

std::cout<<stepCount<<" "<<timeNow<<" "<<VCSum(&vSum) / nMol<<" "
    <<PropEst(totEnergy)<<" "<<PropEstSig(totEnergy)
    <<" "<<PropEst(kinEnergy)<<" "<<PropEstSig(kinEnergy)
    <<" "<<PropEst(pressure)<<" "<<PropEstSig(pressure)<<"\n";
 }; 

double SingleStep(FILE *fp){
    stepCount++;    
    timeNow = stepCount* deltaT;

    clock_t start = clock();


    LeapfrogStep(1);
    ApplyBoundaryCond();
    ComputeForces();
    LeapfrogStep(2);

    EvalProps();
    AccumProps(1);
    for (int i = 0; i < nMol; i++) {
        dataCoords.at(stepCount - 1).push_back({ mol[i].coordinates.x, mol[i].coordinates.y });
    }
    if(stepCount % stepAvg == 0){
        AccumProps(2);
        if(stepCount>=stepEquil && (stepCount - stepEquil) % stepVel == 0){
            EvalVelDist();
        }
        //PrintSummary(fp);
        AccumProps(0);
    }
    clock_t end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    return time_taken;

};


void SetupJob(){
    stepCount = 0;
    countVel = 0;
    InitCoords(); // initial coordinates
    InitVels(); // initial velocities
    InitAccels(); // initial accelerations
    AccumProps(0);
};

void GetNameList(const char* fd, std::vector<KeyValue>* data) {
    std::ifstream file(fd);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << "\n";
        return;
    }

    std::string line;
    const std::string pattern = "initUcell";

    while (std::getline(file, line)) {
        KeyValue kv;

        if (line.substr(0, pattern.size()) == pattern) {
            line.erase(0, pattern.size());
            kv.key = pattern;
            size_t pos = line.find_first_not_of(" \t");
            if (pos != std::string::npos) {
                line.erase(0, pos);
                kv.value = line;
                (*data).push_back(kv);
            }
        } else {
            size_t pos = line.find_first_of(" \t");
            if (pos != std::string::npos) {
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

void PrintNameList (std::vector<KeyValue> *data){
    for(const auto& kv : (*data)) {
        std::cout << kv.key << " " << kv.value << "\n";
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
    hFunction = 0;
    for(j = 0; j<sizeHistVel; j++){
        if(histVel[j] > 0){
            hFunction += histVel[j]*log(histVel[j]/((j+0.5)*deltaV));
        }
    }
}

void PrintVelDist(FILE *fp){
    real vBin;
    int n;

    printf ("vdist (%.3f)\n", timeNow);
    for (n = 0; n < sizeHistVel; n ++) {
        vBin = (n + 0.5) * rangeVel / sizeHistVel;
        fprintf (fp, "%8.3f %8.3f\n", vBin, histVel[n]);
    }
    fprintf (fp, "hfun: %8.3f %8.3f\n", timeNow, hFunction);
}   