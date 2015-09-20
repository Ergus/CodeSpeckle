CXX = icpc
CXXFLAGS = -O2

MKLFLAGS=-I/u/shared/programs/x86_64/mkl/11.1.3/composer_xe_2013_sp1.3.174/mkl/include -mkl=parallel -lpthread

FILE = main.x
LIBS = speckleauxiliary.o specklemod.o inits.o parser.o

all: $(FILE)

debug: CXXFLAGS = -O0 -DDEBUG -g -Wall
debug: $(FILE)

%.x: %.cc $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(MKLFLAGS)

.SUFFIXES: .o .cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(MKLFLAGS) 

.PHONY: clean
clean:
	rm -rf *.o *.x

.PHONY: test
test: $(FILE)
	./$(FILE)

