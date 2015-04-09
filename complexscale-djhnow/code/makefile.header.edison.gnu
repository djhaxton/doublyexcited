F90 = ftn 
CC = cc 

# (djh doesn't know what is absolutely the best)
OPTS = -DUSE_PETSC -C -check all -warn all,nodec,nounused -traceback
#  -check noarg_temp_created
OPTS = -cpp -w -DUSE_PETSC

OPTS = -parallel -ipo -xAVX -dynamic -cpp -w -DUSE_PETSC 

OPTS = -fast -cpp -w -DUSE_PETSC 


#using this for fast
OPTS = -Ofast -cpp -w -fopenmp -DUSE_PETSC

#using this for debug
OPTS = -DUSE_PETSC -fbounds-check -Wall -Wno-unused -dynamic


LIBBLAS=

INCPETSC= -I$(PETSC_DIR)/include
LIBPETSC= -L$(PETSC_DIR)/lib -lcraypetsc_gnu_complex

INCPETSC= 
LIBPETSC= 

INCS=$(INCPETSC) $(INCSLU)
LIBS=$(LIBPETSC) $(LIBSLU) $(LIBMETIS) $(LIBBLAS)

