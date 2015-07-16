include ../Makefile.in

#==============================================================================
# all object files that must be linked into final executable
OBJS = 
OBJS += main.o 
OBJS += gauleg.o
OBJS += bessel.o
OBJS += ellint.o
OBJS += grnfcn.o

#==============================================================================
# rule for compiling everything
all : main

#==============================================================================
# rule for building executable from object files
main : $(OBJS) $(LIBRARIES)
	$(LOADER) $(LOADOPTS) $(OBJS) \
	$(LIBRARIES) -o $@
	./$@

#==============================================================================
# rule for compiling individual source files into object files
.cpp.o :
	$(COMPILE) -o $@ $<

#==============================================================================
# rule to clean up all output files
clean :
	rm -f *.o main
