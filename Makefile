CC = icpc
CFLAGS = -g -Wall -O2

mklflags=-I/u/shared/programs/x86_64/mkl/11.1.3/composer_xe_2013_sp1.3.174/mkl/include -mkl=parallel -lpthread

files = main.x
llibs = speckleauxiliary.o specklemod.o

all: $(files)

%.x: %.cc $(llibs)
	$(CC) $^ -o $@ $(mklflags) 

.SUFFIXES: .o .cc
.cc.o:
	$(CC) -c $< -o $@ $(mklflags) 

.PHONY: clean
clean:
	rm -rf *.o *.x

.PHONY: test
test: $(files)
	./$(files)

