OBJS = NLC_1D_TFIM.o GenHam.o Lanczos_07.o graphs.o #Lattice_16B.cpp
CC = g++
CFLAGS = -O2 -fopenmp -Wall -Wextra --pedantic
#CFLAGS = -O2 -arch x86_64
#LIBS = -lm -framework veclib

1d.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o 1d.out $(LIBS)

NLC_1D_TFIM.o : NLC_1D_TFIM.cpp CPU/GenHam.h CPU/Lanczos_07.h CPU/simparam.h
	$(CC) $(CFLAGS) -c NLC_1D_TFIM.cpp

GenHam.o: CPU/GenHam.cpp CPU/GenHam.h CPU/Lanczos_07.h
	$(CC) $(CFLAGS) -c CPU/GenHam.cpp

#Lattice_16B.o: Lattice_16B.cpp GenHam.h 
#	$(CC) $(CFLAGS) -c Lattice_16B.cpp

Lanczos_07.o: CPU/Lanczos_07.cpp CPU/GenHam.h CPU/Lanczos_07.h
	$(CC) $(CFLAGS) -c CPU/Lanczos_07.cpp

graphs.o: ../Graphs/graphs.cpp ../Graphs/graphs.h
	$(CC) $(CFLAGS) -c ../Graphs/graphs.cpp
clean :
	rm *.o
