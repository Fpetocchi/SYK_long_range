#---------------------------------------------------------------------------------
#   NOTES:
#
#................................................
#          DIRECTORY STRUCTURE
#................................................
BINDIR= exe
SRCDIR= src
OBJDIR= obj
MODDIR= mod
#................................................
#  RULES for precompliler, i.e. F90->f90
#  and to produce object files.
#................................................
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJDIR)/%.o: $(SRCDIR)/%.f
	$(FC77) $(F77FLAGS) -c $< -o $@
$(SRCDIR)/%.f90: $(SRCDIR)/%.F90
	$(FCp) -D$(archlib) $< $@
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) $(LOCINC) -c $< -o $@
#................................................
#                 FLAGS
#................................................
#CC=/opt/local/bin/gcc
#CXX=/opt/local/bin/g++
CC=gcc
CXX=g++
 QPAR=-fopenmp
 LFLAGS =$(OPTI) $(DEBUG)
 CFLAGS =$(OPTI) -I$(SRCDIR) -std=c++0x -I/usr/local/include # -Wall -Wextra -Wimplicit-interface -fmax-errors=1 -g -fcheck=all -fbacktrace
 LIBS = -lstdc++ -lgsl -lgslcblas -lm
 DEBUG=  -Wall -Wextra -Wimplicit-interface -fmax-errors=1 -g -fcheck=all -fbacktrace
#OPTI= -O3 -fpermissive
#................................................
#             Object files
#................................................
obj1= rungekutta.o harmosc.o
obj2= spin2Dxxz.o

# --- add path to obj files --
OBJ1=$(addprefix $(OBJDIR)/,$(obj1))
OBJ2=$(addprefix $(OBJDIR)/,$(obj2))
#================================================
#              RUN PREPROCESSOR
#================================================

#================================================
#          LIST OF PRODUCED BINARIES
#================================================
serial = harmosc.x spin2Dxxz.x
serialbin = $(addprefix $(BINDIR)/,$(serial))

#====================================================
# MAKE TARGETs:
#        clean - remove all objects, modules and
#	         binaries
#        prune - remove backup files created by
#	         emacs editor
#        info  - print basic information about the
#                host
#      serial  - make all serial programs
#====================================================
.PHONY: clean prune purge info
clean:$
	rm -rf $(OBJDIR)/*.o
	rm -rf  $(BINDIR)/*.x
prune:$
	rm -rf $(SRCDIR)/*~
	rm -rf $(SRCDIR)/*#
#	find . -name '*~' -exec rm -rf {} \;
info: $
	@echo MACHTYPE: $(MACHTYPE) HOSTNAME: $(HOSTNAME) mpi: $(mpi) serial: $(serial)
purge:$
	rm -rf OUTPUT/* ERROR/*

all:$ $(serialbin)

#................................................
#       Fortran file dependencies
#................................................


#................................................
#                LINKING
#................................................
$(BINDIR)/harmosc.x:$ $(OBJ1) # Test program for harmonic oscillator
	$(CXX)  -o $@ $^ $(LFLAGS) $(LIBS)
$(BINDIR)/spin2Dxxz.x:$ $(OBJ2) # Test program for 2D spin
	$(CXX)  -o $@ $^ $(LFLAGS) $(LIBS)
