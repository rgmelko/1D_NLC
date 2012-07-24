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
#include <string>

#include "CPU/Lanczos_07.h"
#include "CPU/GenHam.h"
#include "CPU/simparam.h"
#include "graphs.h"

int main(int argc, char** argv){

    int CurrentArg = 1;
    string InputFile;
    string OutputFile = "output_2d.dat";
    while (CurrentArg < argc)
    {
        if (argv[ CurrentArg ] == string("-i") || argv[ CurrentArg ] == string("--input") )
        {
            InputFile = string(argv[ CurrentArg + 1 ]);
        }
        if (argv[ CurrentArg ] == string("-o") || argv[ CurrentArg ] == string("--output"))
        {
            OutputFile = string(argv[ CurrentArg + 1 ]);
        }
        CurrentArg++;
    }

    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    vector <long double> eVec;

    J=prm.JJ_;
    h=prm.hh_;

    vector< graph > fileGraphs; //graph objects
    
    vector<double> WeightHigh;
    
    ReadGraphsFromFile(fileGraphs, InputFile);

    ofstream fout(OutputFile.c_str());
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
        energy = lancz.Diag(HV, 1, 1, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
        WeightHigh.push_back(energy);
        for (int j = 0; j < fileGraphs[ i ].SubgraphList.size(); j++)
	  WeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * WeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];

        cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<"  ";
	//        cout<<" energy "<<setprecision(12)<<energy<<endl;
	//        cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	    RunningSumHigh += fileGraphs[ i ].LatticeConstant * WeightHigh.back();
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
