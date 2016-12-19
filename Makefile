#
#   Makefile for ddCOSMO
#
#RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp 
RunF77 = gfortran
#FFLAGS = -O1 -fbacktrace -g -march=native -fopenmp -llapack -lblas
FFLAGS = -O0 -fbacktrace -fbounds-check -g -march=native -fopenmp -llapack -lblas
#RunF77 = pgfortran
#FFLAGS = -O3 -mp

MODS   = ddcosmo.o
OBJS   = mkrhs.o llgnew.o main.o ddcosmo.o forces.o efld.o iefpcm.o compute_forces.o
#
all:    $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) -o main.exe $(OBJS)
#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(OBJS) *.exe *.mod
