PUB     = /class/public/cse260fa10
#include $(PUB)/Arch/arch.gnu
include arch.gnu

#include ./arch.mac

# If you want to compile for single precision,
# uncomment the next line
# ARITH = -DFLOAT

# This generates output about how the
# compiler vectorized the code
# We  suggest using level 2 (the integer after "verbose=")
# See the gcc manual for the other levels of output: levels 1-7
# http://gcc.gnu.org/onlinedocs/gcc-4.4.2/gcc/Debugging-Options.html#Debugging-Options
#
OPTIMIZATION += -ftree-vectorize -msse2 -ffast-math

C++FLAGS += $(ARITH) -fopenmp
CFLAGS += $(ARITH) -fopenmp -Wall

ifdef debug
CFLAGS += -DDEBUG
endif

#DEBUG += -DDEBUG

#Uncomment the next line  if you want to add symbol table
# information for gdb
#DEBUG += -g


app:		apf

OBJECTS = apf.o solve.o splot.o cmdline.o printMat.o Timer.o

apf:	        $(OBJECTS) 
		$(CLINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

clean:	
	$(RM) *.o apf;
	$(RM) core;
