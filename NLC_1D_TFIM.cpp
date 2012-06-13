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

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat");

    for (int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph

        GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList); 

        LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
        HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
        energy = lancz.Diag(HV, 1, prm.valvec_); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors

        cout<<"graph #"<<i;
        cout<<" energy "<<setprecision(12)<<energy<<endl;

    }

    return 0;

}
