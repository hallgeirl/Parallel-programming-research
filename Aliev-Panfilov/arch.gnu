#########################################################################
#									#
# Sample makefile header for running with Gnu compilers  		#
#  The makefile targets are appended to  the end of this file		#
#	 Don't change anything that comes before the targets 		#
#									#
#									#
#########################################################################

RM		= rm -f
LN		= ln -s
ECHO		= echo


C++ 		= g++
#CC		= gcc
AR		= ar
RANLIB		= ranlib
C++LINK		= $(C++)
CLINK		= $(CC)



ARCH_FLAGS      = -DLINUX 
WARNINGS        = 

C++FLAGS        += $(INCLUDES) $(ARCH_FLAGS) $(WARNINGS) \
                  $(XTRAFLAGS) $(DEBUG)

CFLAGS		+= $(INCLUDES) $(ARCH_FLAGS) $(WARNINGS) \
                  $(XTRAFLAGS) $(DEBUG)

FFLAGS		= $(ARCH_FLAGS) -O2 -fno-second-underscore -ff90 -fugly-complex




ARFLAGS		= ru


LDFLAGS		= $(WARNINGS) $(OPTIMIZATION) $(DEBUG)
LDLIBS		= -lm -lpthread


ARCH_HAS_X	= arch_has_X



#########################################################################
# End of the System dependent prefix
#########################################################################


#########################################################################
#									#
# Suffixes for compiling most normal C++ and  C files		        #
#									#
#########################################################################

.SUFFIXES:
.SUFFIXES: .C .cpp .c .o

.C.o:
		@$(ECHO)
		@$(ECHO) "Compiling Source File --" $<
		@$(ECHO) "---------------------"
		$(CC) $(C++FLAGS) -c $<
		@$(ECHO)

.cpp.o:
		@$(ECHO)
		@$(ECHO) "Compiling Source File --" $<
		@$(ECHO) "---------------------"
		$(C++) $(C++FLAGS) -c $<
		@$(ECHO)



.c.o:
		@$(ECHO)
		@$(ECHO) "Compiling Source File --" $<
		@$(ECHO) "---------------------"
		$(CC) $(CFLAGS) -c $<
		@$(ECHO)

