include ../Makefile.in
#include $(LAPACK_DIR)/make.inc

all: gaussleg

#==============================================================================
# executable
gaussleg: gaussleg.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) gaussleg.o \
	$(LIBRARIES) -o $@
	./$@

#==============================================================================
# dependencies
.cpp.o:
	$(COMPILE) -o $@ $<

#==============================================================================
# clean
clean:
	rm -f *.o gaussleg


