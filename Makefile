#
#   Makefile for ddCOSMO
#
RunF77 = ifort -mkl=parallel
FFLAGS = -O3 -xHost -qopenmp 
#RunF77 = gfortran
#FFLAGS = -O3 -march=native -fopenmp -llapack -lblas

MODS   = ddcosmo.o
OBJS   = mkrhs.o llgnew.o main.o ddcosmo.o forces.o efld.o iefpcm.o compute_forces.o numgrad.o fmm_dummy.o\
	jacobi_diis.o cosmo.o gmres.o matvec.o pcm.o debug.o
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
