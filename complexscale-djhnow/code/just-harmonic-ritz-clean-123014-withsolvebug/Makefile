
# for checkmodules to work, all .o instructions must have object.o:  etc.  not object.o : etc.

F90=gfortran 




FFLAGS= -x f95-cpp-input -ffree-line-length-none -ffast-math -O3 -Wall

FFLAGS= -x f95-cpp-input -ffree-line-length-none -ffast-math -O3 -fbounds-check -Wall



LIBRARIES= -framework Accelerate

SRCS= hr_eigen.o    harmonic_ritz_complex.o  myprogram.o SineDVR.o 

default: myprogram

clean:
	rm *.o *.mod *~  myprogram *genmod*

myprogram:   $(SRCS) 
	$(F90) -o myprogram $(SRCS) $(LIBRARIES)

harmonic_ritz_complex.o: new_harmonic_ritz.F90 shell.h.f90 harmonic_ritz_complex.F90 hr_eigen.o
	$(F90) $(FFLAGS) -c harmonic_ritz_complex.F90

SineDVR.o: SineDVR.F90
	$(F90) $(FFLAGS) -c SineDVR.F90

hr_eigen.o: hr_eigen.f90 
	$(F90) $(FFLAGS) -c hr_eigen.f90 

myprogram.o: myprogram.f90 harmonic_ritz_complex.o
	$(F90) $(FFLAGS) -c myprogram.f90 

arg.o: arg.c 
	$(CC) $(PGFFLAG) -c arg.c

.F90.o: 
	$(F90) $(FFLAGS) -c $<




