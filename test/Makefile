include ../install/lapack-3.5.0/make.inc

all: main

# directories
LAPACKDIR = ../install/lapack-3.5.0
ATLASDIR = ../install/atlas-3.10.2/lib/lib

# libraries (archives)
LAPACKLIB = $(LAPACKDIR)/liblapack.a
LAPACKELIB = $(LAPACKDIR)/liblapacke.a
BLASLIB = $(LAPACKDIR)/librefblas.a
CBLASLIB = $(ATLASDIR)/libcblas.a
ATLASLIB = $(ATLASDIR)/libatlas.a

LIBRARIES= $(LAPACKELIB) $(LAPACKLIB) $(BLASLIB) $(CBLASLIB) $(ATLASLIB)

# executable
main: main.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS)  main.o \
	$(LIBRARIES) -o $@
	./$@

# addition by preyas: similar to his makefile
#main: main.o
#	$(CC) -o $@ main.o

# dependencies
.cpp.o:
	$(CC) -c $(CFLAGS) -I $(LAPACKDIR)/lapacke/include -o $@ $<

# clean
clean:
	rm -f *.o main


