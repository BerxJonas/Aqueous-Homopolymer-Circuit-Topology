CC = gcc
CFLAGS = -O0 -g -Wall -std=c99
LIBS = -lm -lgsl -lgslcblas

all: circuit_topology_analysis.exe

%.exe: %.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.exe *.o core.*
