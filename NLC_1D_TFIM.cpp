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
  

//  GENHAM HV(fileGraphs.at(15).NumberSites,J,h,fileGraphs.at(15).AdjacencyList); 
  //HV.Bonds_16B(); 

J=1;
h=1;
  // ------------------------
  vector < pair < int,int > > Bond;
  Bond.resize(15);
  Bond[0]  = make_pair( 0, 1);
  Bond[1]  = make_pair( 1, 2);
  Bond[2]  = make_pair( 2, 3);
  Bond[3]  = make_pair( 3, 4);
  Bond[4]  = make_pair( 4, 5);
  Bond[5]  = make_pair( 5, 6);
  Bond[6]  = make_pair( 6, 7);
  Bond[7]  = make_pair( 7, 8);
  Bond[8]  = make_pair( 8, 9);
  Bond[9]  = make_pair( 9,10);
  Bond[10] = make_pair(10,11);
  Bond[11] = make_pair(11,12);
  Bond[12] = make_pair(12,13);
  Bond[13] = make_pair(13,14);
  Bond[14] = make_pair(14,15);
//  Bond[15] = make_pair(15, 0);
  // ------------------------

  GENHAM HV(16,J,h,Bond); 

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
