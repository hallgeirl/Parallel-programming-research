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

ifeq ($CC,icc)
else
endif

CFLAGS += $(ARITH) -Wall -O3

ifdef debug
CFLAGS += -DDEBUG -g
endif

app:		apf_pthreads

OBJECTS = apf.o solve.o splot.o cmdline.o printMat.o Timer.o

apf_pthreads:	        $(OBJECTS) 
		$(CLINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

clean:	
	$(RM) *.o apf_pthreads;
	$(RM) core;
