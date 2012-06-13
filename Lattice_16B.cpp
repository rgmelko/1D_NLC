#include "GenHam.h"

// 1D 16 site chain (PBC)
//----------------------------------------------------------
void GENHAM::Bonds_16B(){

  Bond.resize(16);

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
  Bond[15] = make_pair(15, 0);


}//MakeBonds

