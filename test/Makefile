include ../Makefile.in
#include $(LAPACK_DIR)/make.inc

all: main

#==============================================================================
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


