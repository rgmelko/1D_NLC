#include "GPUmagnetization.h"

__host__ void GPUmagnetization( int HowMany, d_hamiltonian* & Ham, parameters* & data,  double** & Eigenvectors, double* & chi)
{
    int VecSize[ HowMany ];
    double** result = (double**)malloc( HowMany * sizeof(double*) ); 
    cudaError_t status[ HowMany ];
    cudaStream_t stream[ HowMany ];

    for( int i = 0; i < HowMany; i++)     
    {

        VecSize[i] = Ham[i].sectorDim;
        
        status[i] = cudaStreamCreate( &stream[i] );
        if( status[i] != cudaSuccess )
        {
            cout<<"Error creating "<<i<<"th stream in magnetization: "<<cudaGetErrorString( status[i] )<<endl;
        }

        status[i] = cudaMalloc( &result[i], VecSize[i]*sizeof(double));
        if( status[i] != cudaSuccess )
        {
            cout<<"Error allocating "<<i<<"th result in magnetization: "<<cudaGetErrorString( status[i] )<<endl;
        }

    }

    for( int i = 0; i < HowMany; i++)
    {
        MagnetSquared<<<VecSize[i]/512, 512, 0, stream[i]>>>(VecSize[i], data[i].nsite, Eigenvectors[i], result[i]);
    }


    for( int i = 0; i < HowMany; i++)
    {
        thrust::device_ptr<double> ReducePtr( result[i] );
        chi[ i ] = thrust::reduce(ReducePtr, ReducePtr + VecSize[i]);
    }
    
    for(int i = 0; i < HowMany; i++)
    {
        cudaFree(result[i]);
        cudaStreamDestroy(stream[i]);
    }
    free( VecSize );
    free( result );
}

__global__ void MagnetSquared( int VecSize, int SpinCount, double* Groundstate, double* MagSquared)
{
    int CurrentKet = threadIdx.x + (blockIdx.x * blockDim.x);
    int UpSpins = 0;

    if( CurrentKet < VecSize )
    {
        for( int CurrentSpin = 0; CurrentSpin < SpinCount; CurrentSpin++)
        {
            UpSpins += ( CurrentKet >> CurrentSpin) & 1;
        }
        MagSquared[ CurrentKet ] = ((2*UpSpins) - SpinCount ) * Groundstate[ CurrentKet ] * Groundstate[ CurrentKet ] * ((2*UpSpins) - SpinCount );
    }
}

