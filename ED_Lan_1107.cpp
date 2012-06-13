/*******************************************************************
Exact Diagonalization program (Lanczos and Complete)
for use with Linked Cluster Expansion
Roger Melko, Ann Kallin, June 2012
baesd on a code from November 2007
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

int main(){

  PARAMS prm;
  double N;
  double J;
  double h; 

  N=prm.NN_;
  J=prm.JJ_;
  h=prm.hh_;


  GENHAM HV(16,J,h); 
  HV.Bonds_16B(); 

  //-------------------------EVERYTHING IN THIS BLOCK OF CODE FOR LANCZOS

  LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)

  HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

  lancz.Diag(HV,prm.Neigen_,prm.valvec_); // second parameter: # of eigenvalues to converge
                      // third parameter: 1 for -values only, 2 for vals AND vectors



  return 0;

}
