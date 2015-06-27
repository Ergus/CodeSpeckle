
CC = gcc
CFLAGS = -g -Wall -O2

files=main.x

all: $(files)

%.x: %.cc
	g++ $^ -o $@

.c.o:
	g++ -c $(CFLAGS) $<
