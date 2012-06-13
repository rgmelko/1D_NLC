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
  double J;
  double J2;
  double Q;
  int Sz;
  
  J=prm.JJ_;
  J2=prm.J2_;
  Q=-prm.QQ_;// SIGN HAS TO BE FLIPPED: NOW SIGN AT INPUT IS SAME AS PAPER 
  Sz=prm.Sz_;

  GENHAM HV(16,J,J2,Q,Sz); 
  HV.Bonds_16B(); 

  //-------------------------EVERYTHING IN THIS BLOCK OF CODE FOR LANCZOS

  LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)

  HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

  lancz.Diag(HV,prm.Neigen_,prm.valvec_); // second parameter: # of eigenvalues to converge
                      // third parameter: 1 for -values only, 2 for vals AND vectors



  return 0;

}
