#include "cuda.h"
#include "iostream"
#include "thrust/device_ptr.h"
#include "thrust/reduce.h"
#include "lanczos.h"

using namespace std;

__host__ void GPUmagnetization( int HowMany, d_hamiltonian* & Ham, parameters* & data, double** & Eigenvectors, double* & chi);

__global__ void MagnetSquared( int VecSize, int SpinCount, double* Groundstate, double* MagSquared);
