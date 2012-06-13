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

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    J=prm.JJ_;
    h=prm.hh_;


    graph graph1;

    vector< graph > fileGraphs;

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat");

    fileGraphs.at(15).print();

    //for( unsigned int currentGraph = 0; currentGraph < fileGraphs.size(); currentGraph++)
    //{
    //    fileGraphs.at(currentGraph).print();
    //    if ( fileGraphs.at(currentGraph).Identifier == 16)
    //    {
    //        graph1 = fileGraphs.at(currentGraph);
    //    }
    //}

    //graph1 = GetGraphFromFile(0,filename);

  //Read in the first line of the graphs file

  GENHAM HV(fileGraphs.at(15).NumberSites,J,h,fileGraphs.at(15).AdjacencyList); 


//  //-------------------------EVERYTHING IN THIS BLOCK OF CODE FOR LANCZOS
//
  LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
//
  HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
//
  lancz.Diag(HV,prm.Neigen_,prm.valvec_); // second parameter: # of eigenvalues to converge
//                      // third parameter: 1 for -values only, 2 for vals AND vectors



  return 0;

}
