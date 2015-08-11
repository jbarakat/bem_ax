include ../Makefile.in

#==============================================================================
# all object files that must be linked into final executable
OBJS = 
OBJS += run.o 
OBJS += gauleg.o gaulog.o
OBJS += bessel.o ellint.o
OBJS += grnfcn.o
OBJS += interp.o

#==============================================================================
# rule for compiling everything
all : run

#==============================================================================
# rule for building executable from object files
run : $(OBJS) $(LIBRARIES)
	$(LOADER) $(LOADOPTS) $(OBJS) \
	$(LIBRARIES) -o $@ -lstdc++
	./$@

#==============================================================================
# rule for compiling individual source files into object files
.cpp.o :
	$(COMPILE) -o $@ $<

#==============================================================================
# rule to clean up all output files
clean :
	rm -f *.o main
