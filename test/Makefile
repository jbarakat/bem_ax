<<<<<<< HEAD
include /home/jmb/Documents/install/lapack-3.5.0/make.inc

all: main

# directories
INSTALLDIR = /home/jmb/Documents/install
LAPACKDIR = $(INSTALLDIR)/lapack-3.5.0
ATLASDIR = $(INSTALLDIR)/atlas-3.10.2/lib/lib

# libraries (archives)
LAPACKLIB = $(LAPACKDIR)/liblapack.a
LAPACKELIB = $(LAPACKDIR)/liblapacke.a
BLASLIB = $(LAPACKDIR)/librefblas.a
CBLASLIB = $(ATLASDIR)/libcblas.a
ATLASLIB = $(ATLASDIR)/libatlas.a

LIBRARIES= $(LAPACKELIB) $(LAPACKLIB) $(BLASLIB) $(CBLASLIB) $(ATLASLIB)

=======
include ../Makefile.in
#include $(LAPACK_DIR)/make.inc

all: main

#==============================================================================
>>>>>>> 75332938999b88e903b9f8b45d84adaa5fc0e707
# executable
main: main.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) main.o \
	$(LIBRARIES) -o $@
	./$@

#==============================================================================
# dependencies
.cpp.o:
	$(COMPILE) -o $@ $<

#==============================================================================
# clean
clean:
	rm -f *.o main


