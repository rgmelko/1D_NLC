/*******************************************************************
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
baesd on a Lanczos code from November 2007
********************************************************************/

#include <utility>
#include <fstream> 
#include <vector> 
#include <math.h>
using namespace std;

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

#include "Lanczos_07.h"
#include "GenHam.h"
#include "lapack.h"
#include "simparam.h"
#include "graphs.h"

int main(){

    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    J=prm.JJ_;
    h=prm.hh_;


    vector< graph > fileGraphs; //graph objects
    
    vector<double> Weight;

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat");

    Weight.push_back(-h); //Weight for site zero
    double RunningSum = Weight[0];
    cout<<RunningSum<<endl;


    for (int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph

        GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList); 

        LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
        HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
        energy = lancz.Diag(HV, 1, prm.valvec_); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors

        Weight.push_back(energy);
        for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++)
            Weight[i] -= fileGraphs.at(i).SubgraphList[j].second * Weight[fileGraphs.at(i).SubgraphList[j].first];

        cout<<"graph #"<<i;
        cout<<" energy "<<setprecision(12)<<energy<<endl;

        cout<<"Weight["<<i<<"] = "<<Weight[i]<<endl;
        RunningSum += Weight[i];
        cout <<"RunningSum = "<< RunningSum;
        cout<<endl;

    }

    return 0;

}
