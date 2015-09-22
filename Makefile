#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...
# Makefile for compile the code speckle
# Master thesis in HPC 2015
# Author Jimmy Aguilar Mena
#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...

SHELL := /bin/bash

CXX = icpc
CXXFLAGS = -O2

# Flags for MKL
MKLFLAGS = -DUMKL -I$(MKLROOT)/include -mkl=parallel -lpthread $(CCFLAGS)

# FLags for Plasma
PLASMAFLAGS = -DUPLASMA -I$(PLASMA_INC) -openmp -DPLASMA_WITH_MKL -I$(HWLOC_HOME)/include $(CCFLAGS)
PLASMALIBS = -L$(PLASMA_LIB) -lplasma -lcoreblasqw -lcoreblas -lquark -llapacke -mkl=parallel -L$(HWLOC_HOME)/lib -lhwloc

# Flags for Magma
MAGMAFLAGS = -DUMAGMA -I$(MAGMA_INC) -I$(CUDADIR)/include -L$(MAGMA_LIB) -L$(CUDADIR)/lib64 -DMAGMA_WITH_MKL -DMAGMA_SETAFFINITY -DADD_ $(CCFLAGS)
MAGMALIBS = -lmagma -lcublas -lcudart  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -openmp

# Check libraries mkl (mandatory), plasma and magma
ifdef MKLROOT
 CXXFLAGS += $(MKLFLAGS)
 LIBS += mkl_solver.o
 ifdef PLASMA_INC	# Now with MKL we check for Plasma, is not mandatory
  ifdef PLASMA_LIB
   CXXFLAGS += $(PLASMAFLAGS)
   CXXFLAGS += $(PLASMALIBS)
   LIBS += plasma_solver.o
  endif
 else
  $(warning No plasma support will be installed,)
  $(warning define PLASMA_INC and PLASMA_LIBS)
 endif
 ifdef MAGMA_INC	       # Now with MKL we check for Plasma, is not mandatory
  ifdef MAGMA_LIB
   CXXFLAGS += $(MAGMAFLAGS)
   CXXFLAGS += $(MAGMALIBS)
   LIBS += magma_solver.o
  endif
 else
  $(warning No Magma support will be installed)
  $(warning define MAGMA_INC and MAGMA_LIBS)
 endif
else
 $(error MKLROOT is undefined, please load MKL)
endif

# Here start the rules.
FILE = main.x
LIBS = auxiliary.o specklemod.o inits.o parser.o

all: $(FILE)

debug: CXXFLAGS = -O0 -DDEBUG -g -Wall
debug: $(FILE)

paralell: load-mpi
parallel: $(FILE)

%.x: %.cc $(LIBS) check-mkl
	$(CXX) $(CXXFLAGS) $^ -o $@

.SUFFIXES: .o .cc
.cc.o: check-mkl
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

# Extra rules
.PHONY: clean test check-mkl load-mpi
clean:
	rm -rf *.o *.x

test: $(FILE)
	./$(FILE) input

load-mpi:
@echo "Lookimg for MPI Compiler"
MPI_RESULT := $(shell which mpicxx 2> /dev/null)
MPI_TEST := $(notdir $(MPI_RESULT))
ifeq ($(MPI_TEST),mpicxx)
  CXX=mpicxx
else
  $(error MPI compiler not found)
endif
