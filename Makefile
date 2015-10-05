#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...
# Makefile for compile the code speckle
# Master thesis in HPC 2015
# Author Jimmy Aguilar Mena
#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...

SHELL := /bin/bash

CXX = g++
CXXFLAGS = -O2 -lpthread -openmp

# Default flags for gcc (is is not working yet)
GCC_MKLFLAGS=-I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

# Default flags for intel
ICC_MKLFLAGS = -I$(MKLROOT)/include -mkl=parallel -lpthread

# Here start the rules.
FILE = speckle.x
FILE_MPI = speckle_mpi.x

LIBS = auxiliary.o histogram.o specklemod.o inits.o parser.o

MPI_RESULT := $(shell which icpc 2> /dev/null)
MPI_TEST := $(notdir $(MPI_RESULT))
ifeq ($(MPI_TEST),icpc)
 CXX=icpc
endif


MPI_RESULT := $(shell which mpicxx 2> /dev/null)
MPI_TEST := $(notdir $(MPI_RESULT))
ifeq ($(MPI_TEST),mpicxx)
 $(info MPI compiler found)
 MPICXX=mpicxx
 MPIFLAGS = -DUMPI
endif

##-----------------------------------------------
# Check libraries mkl (mandatory), plasma and magma
ifdef MKLROOT
# Flags for MKL
MACROS = -DUMKL
ifeq ($(CXX),icpc)
 MKLFLAGS=$(GCC_MKLFLAGS)
else
 MKLFLAGS=$(GCC_MKLFLAGS)
 $(warning No intel compiler available, please load the module intel)
endif
FLAGS = $(MKLFLAGS)
LIBS += mkl_solver.o
endif

ifdef HWLOC_HOME
  MAGMALIBS += -L$(HWLOC_HOME)/lib -lhwloc
endif

ifdef PLASMA_INC	# Now with MKL we check for Plasma, is not mandatory
 ifdef PLASMA_LIB
  ifdef MKLFLAGS
  # FLags for Plasma
  MACROS = -DUPLASMA -DPLASMA_WITH_MKL
  PLASMAFLAGS = -I$(PLASMA_INC) -openmp -I$(HWLOC_HOME)/include $(CCFLAGS)
  PLASMALIBS = -L$(PLASMA_LIB) -lplasma -lcoreblasqw -lcoreblas -lquark
  FLAGS += $(PLASMAFLAGS)
  LIBS += plasma_solver.o
  endif
endif
else
 $(info No plasma support will be installed,)
 $(info define PLASMA_INC and PLASMA_LIBS)
endif

ifdef MAGMA_INC	       # Now with MKL we check for Plasma, is not mandatory
 ifdef MAGMA_LIB
  ifdef MKLFLAGS
  # Flags for Magma
  MACROS += -DUMAGMA -DMAGMA_WITH_MKL -DMAGMA_SETAFFINITY -DADD_
  MAGMAFLAGS = -DUMAGMA -I$(MAGMA_INC) -I$(CUDADIR)/include  $(CCFLAGS)
  MAGMALIBS =  -DUMAGMA -L$(MAGMA_LIB) -L$(CUDADIR)/lib64 -L$(HWLOC_HOME)/lib -lmagma -lcublas -lcudart
  FLAGS += $(MAGMAFLAGS)
  LIBS += magma_solver.o magma_solver_2stage.o
  endif
 endif
else
 $(info No Magma support will be installed)
 $(info define MAGMA_INC and MAGMA_LIBS)
endif
# The macro var contains now the right support to be
# compiles after checking the enviroments
CXXFLAGS += $(MACROS)
#-----------------------------------------------
# to compile everything ones
debug: CXXFLAGS = -O0 -DDEBUG -g -Wall -DUNIX $(MACROS)

debug_mpi: CXXFLAGS = -O0 -DDEBUG -g -Wall -DUNIX $(MACROS)

mpi: CXXFLAGS += $(MPIFLAGS)
mpi: LIBS += slaves_mpi.o

# Show a message if no MKL loaded
ifdef MKLROOT
 all: $(FILE)
 debug: all
 ifeq ($(MPICXX),mpicxx)
  all: mpi
  # mpi debug version
  debug_mpi: mpi
  # mpi normal version
  mpi: $(FILE_MPI)
 else
  all: 
	$(warning Only serial code will be compiled)
	$(warning load openmpi for paralell version)
 endif
else
 all: 
	$(warning Please load MKL)
endif


$(FILE_MPI): main.cc $(LIBS) slaves_mpi.o
	$(MPICXX) $(CXXFLAGS) $(FLAGS) $^ -o $@ $(MAGMALIBS) $(PLASMALIBS)

$(FILE): main.cc $(LIBS) slaves.o
	$(CXX) $(CXXFLAGS) $(FLAGS) $^ -o $@ $(MAGMALIBS) $(PLASMALIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

slaves_mpi.o: slaves.cc
	$(MPICXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

# Extra rules
.PHONY: clean test cleares mklcheck

clean:
	rm -rf *.o *.x

cleanres:
	rm -rf *.out *.log main.x.*,cn*.btr *.dat

test: $(FILE)
	./$(FILE) input

