include ../Makefile.in
#include $(LAPACK_DIR)/make.inc

all: gauleg

#==============================================================================
# executable
gauleg: gauleg.o $(LIBRARIES)
	$(LOADER) $(LOADOPTS) gauleg.o \
	$(LIBRARIES) -o $@
	./$@

#==============================================================================
# dependencies
.cpp.o:
	$(COMPILE) -o $@ $<

#==============================================================================
# clean
clean:
	rm -f *.o gauleg


