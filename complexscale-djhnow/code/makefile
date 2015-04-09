
include makefile.header

default: doublyexcited solvetest

solvetest: solvetest.o 	hr_eigen.o harmonic_ritz_complex.o SineDVR.o utils2.o
	$(F90) $(LINKOPTS) -o solvetest solvetest.o 	hr_eigen.o harmonic_ritz_complex.o SineDVR.o utils2.o $(LIBS)

SineDVR.o: SineDVR.F90
	$(F90) $(FFLAGS) -c SineDVR.F90

doublyexcited: definitions.INC doublyexcited.o \
	input.o structs.o tooth.o \
	utils.o hammod.o initmod.o lanmod.o \
	quad.o \
	solver.o lamod.o \
	petscsolve.o gplmr.o arg.o \
	hr_eigen.o harmonic_ritz_complex.o utils2.o
	$(F90) $(LINKOPTS) -o doublyexcited doublyexcited.o  \
	input.o structs.o tooth.o \
	utils.o hammod.o initmod.o lanmod.o \
	quad.o \
	solver.o lamod.o \
	petscsolve.o gplmr.o arg.o \
	hr_eigen.o harmonic_ritz_complex.o utils2.o \
	$(LIBS)

doublyexcited.o: doublyexcited.F90 structs.o input.o initmod.o lanmod.o solver.o gplmr.o petscsolve.o utils.o definitions.INC
	$(F90) $(OPTS) -c doublyexcited.F90

gplmr.o: gplmr.F90 structs.o lamod.o hammod.o solver.o petscsolve.o harmonic_ritz_complex.o
	$(F90) $(OPTS) $(INCS) -c gplmr.F90

solver.o: littlesolver.F90 hammod.o structs.o definitions.INC
	$(F90) $(OPTS) $(INCS) -c littlesolver.F90 -o solver.o

lamod.o: lamod.f90 
	$(F90) $(OPTS) -c lamod.f90

solvetest.o: solvetest.f90 harmonic_ritz_complex.o hr_eigen.o
	$(F90) $(OPTS) -c solvetest.f90

hr_eigen.o: hr_eigen.f90 
	$(F90) $(OPTS) -c hr_eigen.f90

harmonic_ritz_complex.o: harmonic_ritz_complex.F90  harmonic_ritz_small.F90 hr_eigen.o
	$(F90) $(OPTS) -c harmonic_ritz_complex.F90

utils.o: utils.F90 structs.o input.o lanmod.o
	$(F90) $(OPTS) -c utils.F90

utils2.o: utils2.f90
	$(F90) $(OPTS) -c utils2.f90

lanmod.o: lanmod.F90 structs.o hammod.o definitions.INC
	$(F90) $(OPTS) -c lanmod.F90

initmod.o: initmod.F90 structs.o tooth.o hammod.o definitions.INC
	$(F90) $(OPTS) -c initmod.F90

hammod.o: hammod.F90 structs.o quad.o definitions.INC
	$(F90) $(OPTS) -c hammod.F90

quad.o: quad.F90 structs.o definitions.INC
	$(F90) $(OPTS) -c quad.F90

tooth.o: tooth.f90 definitions.INC
	$(F90) $(OPTS) -c tooth.f90

structs.o: structs.F90 definitions.INC
	$(F90) $(OPTS) -c structs.F90

input.o: input.F90 initmod.o definitions.INC
	$(F90) $(OPTS) -c input.F90

superlusolve.o: superlusolve.f90
	$(F90) $(OPTS) $(INCS) -c superlusolve.f90

petscsolve.o: solver.o petscsolve.F90
	$(F90) $(OPTS) $(INCS) -c petscsolve.F90

arg.o: arg.c
	$(CC) $(COPTS) $(CINCS) -c arg.c

clean:
	rm *.o *.mod doublyexcited solvetest *genmod.f90


