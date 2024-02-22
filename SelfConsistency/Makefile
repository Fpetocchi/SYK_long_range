SHELL := /bin/bash


#
# Default vars------------------------------------------------------------------
# Set ifort as a default if F90 is not set
FC = ifort
# Set compilation of the Group theory lib. false by default
GT = F
# Add a suffix for the main executable. NULL by default
SUFFIX =
ifdef PAD
  SUFFIX = _$(PAD)
endif


#
# Project files-----------------------------------------------------------------
EXECDIR = bin
SRCDIR = src
INCLDIR = include
SRCDIRLIB = src/libs
TSTSDIR = tests
# Main code
SRCS = linalg.f90 fourier_transforms.f90
OBJS = $(SRCS:.f90=.o) main.o
EXE  = SYK


#
# Compiler flags----------------------------------------------------------------
ifeq ($(FC),ifort)
  FFLAGS = -assume byterecl -qopenmp -liomp5 -mkl=sequential -fpp -mcmodel=large -llapack -I$(INCLDIR)
  DBGFFLAGS = -O0 -check bounds -traceback -fpe3 -check all -fp-model precise -fp-stack-check -p -g -warn -debug extended -warn interface --heap-arrays -check all -check noarg_temp_created -qdiag-enable:sc3 -qdiag-enable:sc-parallel3 -C -D_verb
  RELFFLAGS = -ftz -O3
  # -xAVX -xHost  -axCORE-AVX2
endif
ifeq ($(FC),gfortran)
  FFLAGS = -openmp -cpp -ffree-line-length-none -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -llapack -I$(INCLDIR)
  DBGFFLAGS = -O0 -p -g -fbacktrace -Wline-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -fcheck=bounds -Wall  -D_verb #-Wcharacter-truncation
  RELFFLAGS = -cpp -ffree-line-length-none -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -mcmodel=medium
endif


#
# Debug build settings----------------------------------------------------------
DBGDIR = debug
DBGFFLAGS += -module $(DBGDIR)
# Main code
DBGEXE = $(EXECDIR)/$(EXE)_debug$(SUFFIX)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))


#
# Release build settings--------------------------------------------------------
RELDIR = release
RELFFLAGS += -module $(RELDIR)
# Main code
RELEXE = $(EXECDIR)/$(EXE)$(SUFFIX)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))


#
# All rules---------------------------------------------------------------------
.PHONY: all debug release test prepREL prepDBG remake clean


#
# Default build-----------------------------------------------------------------
all: prepREL release


#
# Debug rules-------------------------------------------------------------------
debug: prepDBG $(DBGEXE)

# Modules
$(DBGDIR)/linalg.o: $(SRCDIR)/linalg.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/fourier_transforms.o: $(SRCDIR)/fourier_transforms.f90 $(DBGDIR)/linalg.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

#$(DBGDIR)/bubbles.o: $(SRCDIR)/bubbles.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/input_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o
#	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

# Module contaier, main and post-processing
$(DBGDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) -c $(FFLAGS) -D$(SUFFIX) $(DBGFFLAGS) -o $@ $<

# Executables
$(DBGEXE): $(DBGOBJS)
	$(FC) $(FFLAGS) $(DBGFFLAGS) -o $(DBGEXE) $^ $(DBGLIB)

#
# Release rules-----------------------------------------------------------------
release: prepREL $(RELEXE)

# Modules
$(RELDIR)/linalg.o: $(SRCDIR)/linalg.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/fourier_transforms.o: $(SRCDIR)/fourier_transforms.f90 $(RELDIR)/linalg.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

#$(RELDIR)/bubbles.o: $(SRCDIR)/bubbles.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/input_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o
#	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

# Module contaier, main and post-processing
$(RELDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) -c $(FFLAGS) -D$(SUFFIX) $(RELFFLAGS) -o $@ $<

# Executables
$(RELEXE): $(RELOBJS)
	$(FC) $(FFLAGS) $(RELFFLAGS) -o $(RELEXE) $^ $(RELLIB)

#
# General rules-----------------------------------------------------------------
prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(EXECDIR)
	#@echo $(OBJS)
	#@echo $(DBGOBJS)

prepREL:
	@touch $(SRCDIR)/main.f90
	@paste <( echo "character(len=256),parameter :: revision = \"" ) <( git rev-parse HEAD  ) <( echo "\"") > $(SRCDIR)/revision_SelfCons.inc
	@mkdir -p $(RELDIR_GT) $(RELDIR) $(INCLDIR) $(EXECDIR)

prepDBG:
	@touch $(SRCDIR)/main.f90
	@paste <( echo "character(len=256),parameter :: revision = \"" ) <( git rev-parse HEAD  ) <( echo "\"") > $(SRCDIR)/revision_SelfCons.inc
	@mkdir -p $(DBGDIR_GT) $(DBGDIR) $(INCLDIR) $(EXECDIR)

moveTSTS:
	@mkdir -p $(TSTSDIR)
	@cp $(EXECDIR)/* $(TSTSDIR)

remake: clean all

clean:
	rm -rf $(EXECDIR) $(TSTSDIR) $(INCLDIR) $(DBGDIR) $(RELDIR) $(SRCDIR)/revisionSelfCons.inc
