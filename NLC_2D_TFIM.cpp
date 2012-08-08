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
#include "CPU/magnetization.h"
#include "../Graphs/graphs.h"

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
    double chi;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    vector <long double> eVec;

    J=prm.JJ_;
    h=prm.hh_;

    vector< Graph > fileGraphs; //graph objects
    
    vector<double> EnergyWeightHigh;
    vector<double> MagnetWeightHigh;
    
    ReadGraphsFromFile(fileGraphs, InputFile);

    ofstream fout(OutputFile.c_str());
    fout.precision(10);
    cout.precision(10);
    
    J=1.;
    
    for(double hh = 1; hh < 3; hh += 2.){
      h = hh;
      
      EnergyWeightHigh.push_back(-h); //Weight for site zero
      MagnetWeightHigh.push_back(1.);
      double EnergyRunningSumHigh = EnergyWeightHigh[0];      
      double MagnetRunningSumHigh = MagnetWeightHigh[0];      
      
      for (unsigned int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph
	
	
	//---High-Field---
	    GENHAM HV(fileGraphs.at(i).Order, J, h, fileGraphs.at(i).AdjacencyList, fileGraphs.at(i).LowField); 

        LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
        HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
        energy = lancz.Diag(HV, 1, 2, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
        chi = Magnetization(eVec, fileGraphs.at(i).Order);
        eVec.clear();
        EnergyWeightHigh.push_back(energy);
        MagnetWeightHigh.push_back(chi);

        for (unsigned int j = 0; j < fileGraphs[ i ].SubgraphList.size(); j++)
        {
	        EnergyWeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * EnergyWeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];
	        MagnetWeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * MagnetWeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];
        }
	//        cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
        EnergyRunningSumHigh += fileGraphs[ i ].LatticeConstant * EnergyWeightHigh.back();
	    //cout<<"This cluster's order: "<<fileGraphs[i].Order<<" Lattice Constant: "<<fileGraphs[i].LatticeConstant<<endl;
        //cout<<"Magnetization for this cluster: "<<chi<<endl;
	    //cout<<"Energy for this cluster: "<<energy<<endl;
        //cout<<"Magnetization Weight High: "<<MagnetWeightHigh.back()<<endl;
        if (fabs(fileGraphs[ i ].LatticeConstant * MagnetWeightHigh.back()) > 1.0)
        {
            fout<<"In the "<<fileGraphs.at(i).Identifier<<"th graph, magnetization weight of "<<fileGraphs[ i ].LatticeConstant * MagnetWeightHigh.back()<<endl;
        }
        MagnetRunningSumHigh += fileGraphs[ i ].LatticeConstant * MagnetWeightHigh.back();
	    cout<<"h = "<<hh<<" J = "<<J<<" Energy Running Sum: "<<EnergyRunningSumHigh<<" Magnet Running Sum: "<<MagnetRunningSumHigh<<endl;
      }
      
      fout<<"h= "<<h<<" J= "<<J;	
      fout <<" Energy= "<< EnergyRunningSumHigh<<" Magnet= "<<MagnetRunningSumHigh<<endl;
      
      EnergyWeightHigh.clear();
      MagnetWeightHigh.clear();
      EnergyRunningSumHigh=0;
      MagnetRunningSumHigh=0;
    }
    
    fout.close();    
    return 0;
}
