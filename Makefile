#
#   Makefile for ddCOSMO
#
###RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp 
RunF77 = gfortran
#FFLAGS = -O1 -fbacktrace -g -march=native -fopenmp -llapack -lblas
FFLAGS =  -g -O0 -fbacktrace -fbounds-check -march=native -llapack -lblas -lstdc++
#RunF77 = pgfortran
#FFLAGS = -O3 -mp
#FFLAGS = -O3 -march=native -fopenmp -llapack -lblas

MODS   = ddcosmo.o
OBJS   = mkrhs.o llgnew.o main.o ddcosmo.o forces.o efld.o iefpcm.o compute_forces.o numgrad.o fmm.o\
	jacobi_diis.o cosmo.o gmres.o matvec.o pcm.o debug.o
#
all:    $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) -o main.exe $(OBJS) /home/gatto/RWTH/scalFMM/Build/lib/Release/libscalfmm.a

fmm.o: /home/gatto/RWTH/scalFMM/Tests/Kernels/fmmalone.cpp
	g++ -I/home/gatto/RWTH/scalFMM/Build/Src -I/home/gatto/RWTH/scalFMM/Src -std=c++11 -o fmm.o -c /home/gatto/RWTH/scalFMM/Tests/Kernels/fmmalone.cpp
#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(OBJS) *.exe *.mod
