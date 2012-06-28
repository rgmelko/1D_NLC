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

    Array<l_double,1> eVec;

    J=prm.JJ_;
    h=prm.hh_;


    vector< graph > fileGraphs; //graph objects
    
    vector<double> WeightHigh;

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat");


    ofstream fout("output_1D.dat");
    fout.precision(10);
    cout.precision(10);
    
    J=1;
    
    for(int hh=1; hh<10; hh++){
      h = hh;
      
      WeightHigh.push_back(-h); //Weight for site zero
      double RunningSumHigh = WeightHigh[0];      
      
      for (int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph
	
	
	//---High-Field---
	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,fileGraphs.at(i).LowField); 

        LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
        HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
        energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors

        WeightHigh.push_back(energy);
        for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++)
	  WeightHigh.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightHigh[fileGraphs.at(i).SubgraphList[j].first];

        cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<"  ";
	//        cout<<" energy "<<setprecision(12)<<energy<<endl;
	//        cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	RunningSumHigh += WeightHigh.back();
        cout <<"RunningSumHigh = "<< RunningSumHigh;
        cout<<endl;
	
      }
      
      fout<<"h= "<<h<<" J= "<<J;	
      fout <<" Energy= "<< RunningSumHigh<< endl<<endl;
      
      WeightHigh.clear();
      RunningSumHigh=0;
    }
    
    
    return 0;
    
    fout.close();

}
