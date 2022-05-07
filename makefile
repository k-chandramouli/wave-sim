#-----------------------------------------------------------------------
# This is the makefile for the wave_driver program.  
#
# Possible targets for this descriptor file are wave_exec (the
# default optimized executable), wave_dbg (a debugging version),
# clean (to remove *.o and *.mod files), and realclean to also
# remove the executables.
#-----------------------------------------------------------------------

# The following are considered macros.  We'll use them
# to specify a shell and the compiler and loader programs.

SHELL = /bin/sh
FCOMP = gfortran -fconvert='big-endian'
FLDR = gfortran -fconvert='big-endian'

# ar and ranlib commands

AR = ar
RANLIB = ranlib

# Compiler and loader-option flags are typically listed in
# separate macros.  The macro with -O3 specifies optimization,
# while the macro with -g -C is used for debugging and
# runtime checks of array bounds.


FFLAGS_OPT = -O3
FFLAGS_DBG = -g -fbounds-check

LDRFLAGS = 

# Library paths and libraries can be listed here.

LIBDIR =  -L/usr/local/lib
# LIBS = -llapack -lblas
LIBS = -lopenblas -lpthread

#-----------------------------------------------------------------------

# The following macro lists all of the object files that are
# needed to build the executable.  The "\" signifies
# that the line is continued.

OBJS = wave_kind_mod.o wave_inp_mod.o wave_io_mod.o wave_data_mod.o \
	   wave_timer.o wave_fdlap_mat.o wave_driver.o

MODS = wave_kind_mod.mod wave_inp_mod.mod wave_io_mod.mod wave_data_mod.mod

#-----------------------------------------------------------------------

# The first dependency list is the default, so it should be
# used for the final executable.  Other dependencies may get
# invoked in the process of creating the final executable.
# Note that $(macro) evaluates a macro or an environment variable,
# and the list of commands or rules follow the dependency.

wave_exec : FFLAGS = $(FFLAGS_OPT)
wave_exec : $(OBJS)
	@echo "Creating "$@" in directory "$(PWD)"."
	$(FLDR) -o $@ $(LDRFLAGS) $(OBJS) $(LIBDIR) $(LIBS)

# The following dependency is similar, but it conditionally
# replaces the value of FFLAGS with FFLAGS_DBG when
# processing the rules for the target wave_dbg.

wave_dbg : FFLAGS = $(FFLAGS_DBG)
wave_dbg : $(OBJS)
	@echo "Creating "$@" in directory "$(PWD)"."
	$(FLDR) -o $@ $(LDRFLAGS) $(OBJS) $(LIBDIR) $(LIBS)

#-----------------------------------------------------------------------

# The following are static pattern rules that provide
# instructions for creating object files from the fortran
# source files.

$(OBJS) : %.o : %.f90
	$(FCOMP) $(FFLAGS) -c $<

$(MODS) : %.mod : %.f90
	$(FCOMP) $(FFLAGS) -c $<

# We also need to list the module dependencies.

wave_io_mod.o: wave_kind_mod.mod
wave_inp_mod.o: wave_kind_mod.mod
wave_data_mod.o: wave_kind_mod.mod wave_inp_mod.mod
wave_fdlap+mat.o: wave_kind_mod.mod wave_inp_mod.mod wave_data_mod.mod
wave_time.o: wave_kind_mod.mod
wave_driver.o: wave_kind_mod.mod wave_io_mod.mod wave_inp_mod.mod \
          wave_data_mod.mod

#-----------------------------------------------------------------------

# Specifying 'clean' dependencies is also handy.  With no
# listed dependencies, clean is considered a 'phony' target.
# The dash in front of a command means to continue even if that
# command returns a nonzero error code.

clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)"."

realclean: clean
	-rm wave_exec wave_dbg
	@echo "Cleaned executables from "$(PWD)"."

