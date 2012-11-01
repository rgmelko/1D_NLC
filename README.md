NLC
======

A Numerical Linked Cluster Expansion program for the transverse field Ising model in one and two dimensions.

Developed by [Roger Melko](http://github.com/rgmelko), [Ann Kallin](http://github.com/akallin), and [Katie Hyatt](http://github.com/kshyatt). 

This project aims to provide a free implementation of the [Numerical Linked Cluster Expansion](http://arxiv.org/abs/0706.3254) technique developed by Rajiv Singh, Tyler Bryant, and Marcos Rigol. It is written in C++. The application has two possible run-modes: CPU only, and CPU-GPU. In the second, large clusters are handled using a GPGPU code written in [CUDA](http://www.nvidia.com/object/cuda_home_new.html). 

Documentation can be found [here](http://kshyatt.github.com/1D_NLC/).

Installation
============

You will need to install a few libraries to run the NLC code. If you would like to use the GPU features, you should install the [CUDA toolkit](http://www.nvidia.com/content/cuda/cuda-downloads.html), following the instructions provided by NVIDIA for your operating system and architecture. You should also install the LAPACK and BLAS libraries, which can be found in both apt and yum. If you are using Mac OS/X, you will need to modify the makefiles. Remove `-lblas` and `-llapack` from the `LIBS` line, and add instead `-framework Accelerate`. You will also need to add `-arch=x86_64` to `CFLAGS`. If you want to use the GPU code, you should download it from here: b369fb264fe2de074c04deef9afd69fce023580d . You may need to change the paths in `NLC_*_TFIM.cu` and in `makefile2d` to use this code. Simply change the `../CUDA/Lanczos` to the folder where you downloaded `hamiltonian.*`, `lanczos.*`, etc. 

Running NLC
===========

NLC takes several arguments from `param.dat`, including whether you would like to find the groundstate of the Hamiltonian (not necessary if you only care about energies), the Sz spin sector you are interested in (not relevant for the TF Ising model, as it is not Sz preserving), and the value of `J`, the coupling constant for the bond operators. NLC iterates over several values of `h`, the magnetic field strength. You can easily modify these in `NLC_*_TFIM.*`. NLC reads in the clusters it will operate on. You can specify the file which contains these clusters in the source file for NLC. We plan to allow you to pass the filename from the command line in a future release. If you would like to turn GPU usage on, simply pass `-g` or `--gpu` as a command line argument to your executable. GPU selection should be available in a future release. NLC will dump the final energies it calculates into a file whose name you can specify in `NLC_*_TFIM.*`. If you want to use/create your own clusters, you can use the [graphs](http://github.com/rgmelko/Graphs) utility. Presently this only supports graphs on a square lattice, but further releases should extend this capability.

TODO
====

- 2D TFIM: should we have a different funcion for Lanczos and Heisenberg?
- Change Lanczos convergence to depend on groundstate
- Extend the Hamiltonian code to other models (Heisenberg, J1-J2, XY, etc)
- Allow the user to specify which GPU(s) they would like to use when they launch the program
- Let the user compile with MPI 

------------------------------------------------------------
The original Lanczos was branched from this SVN repository:

URL: svn://saskeram.cmt.uwaterloo.ca/Users/repos/projects/JQ_lanczos/trunk
Repository Root: svn://saskeram.cmt.uwaterloo.ca/Users/repos/projects/JQ_lanczos
Repository UUID: af8637f1-0a1e-4266-a3bc-430f206c0bf5
Revision: 51
Node Kind: directory
Schedule: normal
Last Changed Author: rgmelko
Last Changed Rev: 45
Last Changed Date: 2009-10-16 07:01:23 -0400 (Fri, 16 Oct 2009)
