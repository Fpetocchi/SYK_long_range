#
# Compiler flags
#
#FC     = gfortran
#FFLAGS = -openmp -cpp -ffree-line-length-none -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -llapack
#DBGFFLAGS = -O0 -p -g -fbacktrace -Wline-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -fcheck=bounds -Wall #-Wcharacter-truncation
#RELFFLAGS = -cpp -ffree-line-length-none -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -mcmodel=medium
#
FC     = ifort
FFLAGS = -assume byterecl -qopenmp -liomp5 -mkl=sequential
DBGFFLAGS = -O0 -check bounds -traceback -fpe0 -check all -fp-model precise -fp-stack-check -p -g -warn -debug extended -warn interface
RELFFLAGS = -xHost -O3 -ftz -mcmodel=large


#
# Project files-----------------------------------------------------------------
EXECDIR = bin
SRCDIR = src
SRCS = linalg.f90 parameters.f90 global_vars.f90 utils_misc.f90 utils_fields.f90 file_io.f90 interactions.f90 crystal.f90 fourier_transforms.f90 bubbles.f90 self_energy.f90
OBJS = $(SRCS:.f90=.o)
EXE  = SelfConsistency


#
# Debug build settings----------------------------------------------------------
DBGDIR = debug
DBGEXE = $(EXECDIR)/$(EXE)_debug
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGFFLAGS += -module $(DBGDIR)


#
# Release build settings--------------------------------------------------------
RELDIR = release
RELEXE = $(EXECDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELFFLAGS += -module $(RELDIR)


#
# All rules---------------------------------------------------------------------
.PHONY: all clean debug prep release remake


#
# Default build-----------------------------------------------------------------
all: prep release


#
# Debug rules-------------------------------------------------------------------
debug: prep $(DBGEXE)

#$(DBGEXE): $(DBGOBJS)
#	$(FC) $(FFLAGS) $(DBGFFLAGS) -o $(DBGEXE) $^
#
# Modules
$(DBGDIR)/linalg.o: $(SRCDIR)/linalg.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/utils_misc.o: $(SRCDIR)/utils_misc.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/crystal.o: $(SRCDIR)/crystal.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/fourier_transforms.o: $(SRCDIR)/fourier_transforms.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/parameters.o: $(SRCDIR)/parameters.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/global_vars.o: $(SRCDIR)/global_vars.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/utils_fields.o: $(SRCDIR)/utils_fields.f90 $(DBGDIR)/parameters.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/file_io.o: $(SRCDIR)/file_io.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/parameters.o $(DBGDIR)/utils_fields.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/interactions.o: $(SRCDIR)/interactions.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/parameters.o $(DBGDIR)/global_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/bubbles.o: $(SRCDIR)/bubbles.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/global_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/self_energy.o: $(SRCDIR)/self_energy.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/global_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

#
# Module contaier and main
$(DBGDIR)/module_container.o: $(SRCDIR)/module_container.f90 $(DBGOBJS)
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/main.o: $(SRCDIR)/main.f90 $(DBGDIR)/module_container.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGEXE): $(DBGDIR)/main.o $(DBGDIR)/module_container.o
	$(FC) $(FFLAGS) $(DBGFFLAGS) -o $(DBGEXE) $^


#
# Release rules-----------------------------------------------------------------
release: prep $(RELEXE)

#$(RELEXE): $(RELOBJS)
#	$(FC) $(FFLAGS) $(RELFFLAGS) -o $(RELEXE) $^
#
# Modules
$(RELDIR)/linalg.o: $(SRCDIR)/linalg.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/utils_misc.o: $(SRCDIR)/utils_misc.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/crystal.o: $(SRCDIR)/crystal.f90 $(RELDIR)/utils_misc.o $(RELDIR)/linalg.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/fourier_transforms.o: $(SRCDIR)/fourier_transforms.f90 $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/parameters.o: $(SRCDIR)/parameters.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/global_vars.o: $(SRCDIR)/global_vars.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/utils_fields.o: $(SRCDIR)/utils_fields.f90 $(RELDIR)/parameters.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/file_io.o: $(SRCDIR)/file_io.f90 $(RELDIR)/utils_misc.o $(RELDIR)/parameters.o $(RELDIR)/utils_fields.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/interactions.o: $(SRCDIR)/interactions.f90 $(RELDIR)/utils_misc.o $(RELDIR)/parameters.o $(RELDIR)/global_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/bubbles.o: $(SRCDIR)/bubbles.f90 $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/global_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/self_energy.o: $(SRCDIR)/self_energy.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/global_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

#
# Module contaier and main
$(RELDIR)/module_container.o: $(SRCDIR)/module_container.f90 $(RELOBJS)
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/main.o: $(SRCDIR)/main.f90 $(RELDIR)/module_container.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELEXE): $(RELDIR)/main.o $(RELDIR)/module_container.o
	$(FC) $(FFLAGS) $(RELFFLAGS) -o $(RELEXE) $^


#
# Other rules-------------------------------------------------------------------
prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(EXECDIR)

remake: clean all

clean:
	rm -rf $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS) $(DBGDIR) $(RELDIR) $(EXECDIR)
