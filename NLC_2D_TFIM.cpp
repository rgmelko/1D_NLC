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
#include <mpi.h>

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
    double J = 1.;

    vector <long double> eVec;

    vector< graph > fileGraphs; //graph objects
    
    vector<double> WeightHigh;
    
    ReadGraphsFromFile(fileGraphs, InputFile);
    
    if ( rank == 0 ) 
    {
        ofstream fout(OutputFile.c_str());
        fout.precision(10);
        double* Results = (double*)malloc((size - 1)*sizeof(double));
        for( int i = 0; i < size - 1; i++ )
        {
            MPI_Recv(Results + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &status); //grab results from other processes
            fout<<"h = "<<(i + 1)/J<<" J = "<<J<<" Energy = "<<Results[i]<<endl;
        }
        fout.close();
    }
    
    double h = (double)rank/J;
      
    WeightHigh.push_back(-h); //Weight for site zero
    double RunningSumHigh = WeightHigh[0];      
    
    if ( rank )
    {
        for (int i = 1; i < fileGraphs.size(); i++)
        { //skip the zeroth graph
	
	    //---High-Field---
	        GENHAM HV(fileGraphs.at(i).NumberSites, J, h, fileGraphs.at(i).AdjacencyList, fileGraphs.at(i).LowField); 

            LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
    
            HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
            energy = lancz.Diag(HV, 1, 1, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
            WeightHigh.push_back(energy);
    
            for (int j = 0; j < fileGraphs[ i ].SubgraphList.size(); j++)
            {
	            WeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * WeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];
            }

            RunningSumHigh += fileGraphs[ i ].LatticeConstant * WeightHigh.back();
	
        }
        MPI_Send( &RunningSumHigh, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //send resulting sum to controller process

    }

    return 0;

}
