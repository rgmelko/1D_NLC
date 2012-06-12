#include "GenHam.h"

// Square (4x4) conventional 16 site lattice
//----------------------------------------------------------
void GENHAM::Bonds_16B(){

  Bond.resize(Nsite,3);
  //horizontal    vertical
  Bond = 0, 1, 4,   //   This is all you need for the Heisenberg model
         1, 2, 5,   //   Each row indexes two sites associated with ll
         2, 3, 6,   //   
         3, 4, 7,   //     4
         4, 5, 8,   //     |           e.g. site 0 related to site 1 and 4
         5, 6, 9,   //   ( 0 ) - 1
         6, 7, 10,
         7, 8, 11,
         8, 9, 12,
         9, 10, 13,
         10, 11, 14,
         11, 12, 15,
         12, 13, 0,
         13, 14, 1,
         14, 15, 2,
         15, 0, 3;


}//MakeBonds

