#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...
# Makefile for compile the code speckle
# Master thesis in HPC 2015
# Author Jimmy Aguilar Mena
#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...

SHELL := /bin/bash

CXX = icpc
CXXFLAGS = -O2

# Here start the rules.
FILE = main.x
LIBS = auxiliary.o specklemod.o inits.o parser.o

##-----------------------------------------------
# Check libraries mkl (mandatory), plasma and magma


ifdef MKLROOT
# Flags for MKL
 MACROS = -DUMKL
 MKLFLAGS = -I$(MKLROOT)/include -mkl=parallel -lpthread
 FLAGS = $(MKLFLAGS)
 LIBS += mkl_solver.o
else
 $(error MKLROOT is undefined, please load MKL)
endif

ifdef PLASMA_INC	# Now with MKL we check for Plasma, is not mandatory
 ifdef PLASMA_LIB
  # FLags for Plasma
  MACROS = -DUPLASMA -DPLASMA_WITH_MKL
  PLASMAFLAGS = -I$(PLASMA_INC) -openmp -I$(HWLOC_HOME)/include $(CCFLAGS)
  PLASMALIBS = -L$(PLASMA_LIB) -lplasma -lcoreblasqw -lcoreblas -lquark -llapacke -mkl=parallel
  FLAGS += $(PLASMAFLAGS)
  LIBS += plasma_solver.o
endif
else
 $(warning No plasma support will be installed,)
 $(warning define PLASMA_INC and PLASMA_LIBS)
endif

ifdef MAGMA_INC	       # Now with MKL we check for Plasma, is not mandatory
 ifdef MAGMA_LIB
  # Flags for Magma
  MACROS += -DUMAGMA -DMAGMA_WITH_MKL -DMAGMA_SETAFFINITY -DADD_
  MAGMAFLAGS = -DUMAGMA -I$(MAGMA_INC) -I$(CUDADIR)/include  $(CCFLAGS)
  MAGMALIBS =  -DUMAGMA -L$(MAGMA_LIB) -L$(CUDADIR)/lib64 -lmagma -lcublas -lcudart  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -openmp
  FLAGS += $(MAGMAFLAGS)
  LIBS += magma_solver.o magma_solver_2stage.o
 endif
else
 $(warning No Magma support will be installed)
 $(warning define MAGMA_INC and MAGMA_LIBS)
endif

CXXFLAGS += $(MACROS)
##-----------------------------------------------

all: $(FILE)

debug: CXXFLAGS = -O0 -DDEBUG -g -Wall -DUNIX $(MACROS)
debug: $(FILE)

%.x: %.cc $(LIBS)
	$(CXX) $(CXXFLAGS) $(FLAGS) $^ -o $@ $(MAGMALIBS) $(PLASMALIBS)

plasma_solver.o: plasma_solver.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

magma_solver.o: magma_solver.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

magma_solver_2stage.o: magma_solver_2stage.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

mkl_solver.o: mkl_solver.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

specklemod.o: specklemod.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -c $^ -o $@

# Extra rules
.PHONY: clean test load-mpi check-libs check-mkl check-plasma check-magma

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


