OBJS = NLC_1D_TFIM.cpp GenHam.o  Lanczos_07.o lapack.o graphs.o #Lattice_16B.cpp
CC = g++
#CFLAGS = -O2 
CFLAGS = -O2 -arch x86_64
#LIBS = -lm -framework veclib
LIBS = -framework Accelerate

a.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o a.out $(LIBS)

NLC_1D_TFIM.o : NLC_1D_TFIM.cpp GenHam.h Lanczos_07.h lapack.h simparam.h
	$(CC) $(CFLAGS) -c ED_Lan_1107.cpp

GenHam.o: GenHam.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c GenHam.cpp

#Lattice_16B.o: Lattice_16B.cpp GenHam.h 
#	$(CC) $(CFLAGS) -c Lattice_16B.cpp

Lanczos_07.o: Lanczos_07.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c Lanczos_07.cpp

lapack.o: lapack.cpp lapack.h 
	$(CC) $(CFLAGS) -c lapack.cpp

graphs.o: graphs.cpp graphs.h
	$(CC) $(CFLAGS) -c graphs.cpp
clean :
	rm *.o
