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

#include "CPU/Lanczos_07.h"
#include "CPU/GenHam.h"
#include "CPU/simparam.h"
#include "../Graphs/graphs.h"

int main(){

    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    vector<long double> eVec;

    J=prm.JJ_;
    h=prm.hh_;


    vector< Graph* > fileGraphs; //graph objects
    
    vector<double> WeightHigh;

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat", 0);


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
	    GENHAM HV(fileGraphs.at(i)->Order,J,h,fileGraphs.at(i)->AdjacencyList, 0); 

        LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
        HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
        HV.printg();
        energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
        cout<<"Energy: "<<energy<<endl;
        WeightHigh.push_back(energy);
        for (int j = 0; j<fileGraphs.at(i)->SubgraphList.size(); j++)
	  WeightHigh.back() -= fileGraphs.at(i)->SubgraphList[j].second * WeightHigh[fileGraphs.at(i)->SubgraphList[j].first];

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
