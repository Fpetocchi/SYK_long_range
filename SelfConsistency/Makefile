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
SRCS = linalg.f90 parameters.f90 input_vars.f90 utils_misc.f90 fit.f90 utils_fields.f90 file_io.f90 interactions.f90 crystal.f90 \
       fourier_transforms.f90 bubbles.f90 post_processing.f90 self_energy.f90 greens_function.f90 utils_main.f90 module_container.f90
OBJS = $(SRCS:.f90=.o) minimize.o gap_equation.o main.o
EXE  = SelfConsistency
# Post-processing
OBJS_PP1 = $(SRCS:.f90=.o) minimize.o gap_equation.o pp_Akf_builder.o
EXE_PP1  = Akf_builder
OBJS_PP2 = linalg.o utils_misc.o pp_Gaux2Sigma.o
EXE_PP2  = Gaux2Sigma
#OBJS_PP3 = $(SRCS:.f90=.o) minimize.o gap_equation.o pp_dispersion_fit.o
#EXE_PP3  = Hk_fit


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
# Post-processing
DBGEXE_PP1 = $(EXECDIR)/$(EXE_PP1)_debug
DBGOBJS_PP1 = $(addprefix $(DBGDIR)/, $(OBJS_PP1))
DBGEXE_PP2 = $(EXECDIR)/$(EXE_PP2)_debug
DBGOBJS_PP2 = $(addprefix $(DBGDIR)/, $(OBJS_PP2))
# Custom libraries
DBGDIR_GT = debug/libs
ifeq ($(GT),T)
  KLIB = -D_Klib
  DBGLIB = $(INCLDIR)/libkptgen_debug.a $(INCLDIR)/libutils_debug.a $(INCLDIR)/libsym_debug.a $(INCLDIR)/librational_debug.a
endif


#
# Release build settings--------------------------------------------------------
RELDIR = release
RELFFLAGS += -module $(RELDIR)
# Main code
RELEXE = $(EXECDIR)/$(EXE)$(SUFFIX)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
# Post-processing
RELEXE_PP1 = $(EXECDIR)/$(EXE_PP1)
RELOBJS_PP1 = $(addprefix $(RELDIR)/, $(OBJS_PP1))
RELEXE_PP2 = $(EXECDIR)/$(EXE_PP2)
RELOBJS_PP2 = $(addprefix $(RELDIR)/, $(OBJS_PP2))
# Custom libraries
RELDIR_GT = release/libs
ifeq ($(GT),T)
  KLIB = -D_Klib
  RELLIB = $(INCLDIR)/libkptgen.a $(INCLDIR)/libutils.a $(INCLDIR)/libsym.a $(INCLDIR)/librational.a
endif


#
# All rules---------------------------------------------------------------------
.PHONY: all debug release test prepREL prepDBG remake clean


#
# Default build-----------------------------------------------------------------
all: prepREL release


#
# Debug rules-------------------------------------------------------------------
debug: prepDBG DBGGTlib $(DBGEXE) $(DBGEXE_PP1) $(DBGEXE_PP2)

# Modules
$(DBGDIR)/linalg.o: $(SRCDIR)/linalg.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/utils_misc.o: $(SRCDIR)/utils_misc.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/crystal.o: $(SRCDIR)/crystal.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o
	$(FC) -c $(FFLAGS) $(KLIB) $(DBGFFLAGS) -o $@ $< $(DBGKLIB)

$(DBGDIR)/fourier_transforms.o: $(SRCDIR)/fourier_transforms.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/parameters.o: $(SRCDIR)/parameters.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/input_vars.o: $(SRCDIR)/input_vars.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/parameters.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/utils_fields.o: $(SRCDIR)/utils_fields.f90 $(DBGDIR)/parameters.o $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/file_io.o: $(SRCDIR)/file_io.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/parameters.o $(DBGDIR)/utils_fields.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/minimize.o: $(SRCDIR)/minimize.f
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/fit.o: $(SRCDIR)/fit.f90 $(DBGDIR)/minimize.o $(DBGDIR)/utils_misc.o $(DBGDIR)/parameters.o $(DBGDIR)/file_io.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/interactions.o: $(SRCDIR)/interactions.f90 $(DBGDIR)/utils_misc.o $(DBGDIR)/parameters.o $(DBGDIR)/input_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/crystal.o $(DBGDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/bubbles.o: $(SRCDIR)/bubbles.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/input_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/greens_function.o: $(SRCDIR)/greens_function.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/input_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/self_energy.o: $(SRCDIR)/self_energy.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/input_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o $(DBGDIR)/greens_function.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/gap_equation.o: $(SRCDIR)/gap_equation.f90 $(DBGDIR)/parameters.o $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/post_processing.o: $(SRCDIR)/post_processing.f90 $(DBGDIR)/linalg.o $(DBGDIR)/utils_misc.o $(DBGDIR)/crystal.o $(DBGDIR)/parameters.o $(DBGDIR)/input_vars.o $(DBGDIR)/utils_fields.o $(DBGDIR)/file_io.o $(DBGDIR)/fourier_transforms.o $(DBGDIR)/greens_function.o $(DBGDIR)/gap_equation.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

# Module contaier, main and post-processing
$(DBGDIR)/module_container.o: $(SRCDIR)/module_container.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/utils_main.o: $(SRCDIR)/utils_main.f90 $(DBGDIR)/parameters.o $(DBGDIR)/module_container.o
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) -c $(FFLAGS) -D$(SUFFIX) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/pp_Akf_builder.o: $(SRCDIR)/pp_Akf_builder.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

$(DBGDIR)/pp_Gaux2Sigma.o: $(SRCDIR)/pp_Gaux2Sigma.f90
	$(FC) -c $(FFLAGS) $(DBGFFLAGS) -o $@ $<

# Executables
$(DBGEXE): $(DBGOBJS)
	$(FC) $(FFLAGS) $(DBGFFLAGS) -o $(DBGEXE) $^ $(DBGLIB)

$(DBGEXE_PP1): $(DBGOBJS_PP1)
	$(FC) $(FFLAGS) $(DBGFFLAGS) -o $(DBGEXE_PP1) $^ $(DBGLIB)

$(DBGEXE_PP2): $(DBGOBJS_PP2)
	$(FC) $(FFLAGS) $(DBGFFLAGS) -o $(DBGEXE_PP2) $^ $(DBGLIB)

#
# Release rules-----------------------------------------------------------------
release: prepREL RELGTlib $(RELEXE) $(RELEXE_PP1) $(RELEXE_PP2)

# Modules
$(RELDIR)/linalg.o: $(SRCDIR)/linalg.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/utils_misc.o: $(SRCDIR)/utils_misc.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/crystal.o: $(SRCDIR)/crystal.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $< $(KLIB)

$(RELDIR)/fourier_transforms.o: $(SRCDIR)/fourier_transforms.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/parameters.o: $(SRCDIR)/parameters.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/input_vars.o: $(SRCDIR)/input_vars.f90 $(RELDIR)/utils_misc.o $(RELDIR)/parameters.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/utils_fields.o: $(SRCDIR)/utils_fields.f90 $(RELDIR)/parameters.o $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/file_io.o: $(SRCDIR)/file_io.f90 $(RELDIR)/utils_misc.o $(RELDIR)/parameters.o $(RELDIR)/utils_fields.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/minimize.o: $(SRCDIR)/minimize.f
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/fit.o: $(SRCDIR)/fit.f90 $(RELDIR)/minimize.o $(RELDIR)/utils_misc.o $(RELDIR)/parameters.o $(RELDIR)/file_io.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/interactions.o: $(SRCDIR)/interactions.f90 $(RELDIR)/utils_misc.o $(RELDIR)/parameters.o $(RELDIR)/input_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/crystal.o $(RELDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/bubbles.o: $(SRCDIR)/bubbles.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/input_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/greens_function.o: $(SRCDIR)/greens_function.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/input_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/self_energy.o: $(SRCDIR)/self_energy.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/input_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o $(RELDIR)/greens_function.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/gap_equation.o: $(SRCDIR)/gap_equation.f90 $(RELDIR)/parameters.o $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/post_processing.o: $(SRCDIR)/post_processing.f90 $(RELDIR)/linalg.o $(RELDIR)/utils_misc.o $(RELDIR)/crystal.o $(RELDIR)/parameters.o $(RELDIR)/input_vars.o $(RELDIR)/utils_fields.o $(RELDIR)/file_io.o $(RELDIR)/fourier_transforms.o $(RELDIR)/greens_function.o $(RELDIR)/gap_equation.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

# Module contaier, main and post-processing
$(RELDIR)/module_container.o: $(SRCDIR)/module_container.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/utils_main.o: $(SRCDIR)/utils_main.f90 $(RELDIR)/parameters.o $(RELDIR)/module_container.o
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) -c $(FFLAGS) -D$(SUFFIX) $(RELFFLAGS) -o $@ $<

$(RELDIR)/pp_Akf_builder.o: $(SRCDIR)/pp_Akf_builder.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

$(RELDIR)/pp_Gaux2Sigma.o: $(SRCDIR)/pp_Gaux2Sigma.f90
	$(FC) -c $(FFLAGS) $(RELFFLAGS) -o $@ $<

# Executables
$(RELEXE): $(RELOBJS)
	$(FC) $(FFLAGS) $(RELFFLAGS) -o $(RELEXE) $^ $(RELLIB)

$(RELEXE_PP1): $(RELOBJS_PP1)
	$(FC) $(FFLAGS) $(RELFFLAGS) -o $(RELEXE_PP1) $^ $(RELLIB)

$(RELEXE_PP2): $(RELOBJS_PP2)
	$(FC) $(FFLAGS) $(RELFFLAGS) -o $(RELEXE_PP2) $^ $(RELLIB)

#
# GTlib debug rules-------------------------------------------------------------
DBGGTlib:
ifeq ($(GT),T)
	#
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/num_types.o $(SRCDIRLIB)/symlib/num_types.f90
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/utilities.o $(SRCDIRLIB)/symlib/utilities.f90
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/numerical_utilities.o $(SRCDIRLIB)/symlib/numerical_utilities.f90
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/vector_matrix_utilities.o $(SRCDIRLIB)/symlib/vector_matrix_utilities.f90
	ar ru $(INCLDIR)/libutils_debug.a $(DBGDIR_GT)/num_types.o $(DBGDIR_GT)/utilities.o $(DBGDIR_GT)/numerical_utilities.o $(DBGDIR_GT)/vector_matrix_utilities.o
	ranlib $(INCLDIR)/libutils_debug.a
	#
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/symmetry.o $(SRCDIRLIB)/symlib/symmetry.f90
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/group_theory.o $(SRCDIRLIB)/symlib/group_theory.f90
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/itertools.o $(SRCDIRLIB)/symlib/itertools.f90
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/classes.o $(SRCDIRLIB)/symlib/classes.f90
	ar ru $(INCLDIR)/libsym_debug.a $(DBGDIR_GT)/symmetry.o $(DBGDIR_GT)/group_theory.o $(DBGDIR_GT)/itertools.o $(DBGDIR_GT)/classes.o
	ranlib $(INCLDIR)/libsym_debug.a
	#
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/compare_structures.o $(SRCDIRLIB)/symlib/compare_structures.f90
	ar ru $(INCLDIR)/libcomparestructs_debug.a $(DBGDIR_GT)/compare_structures.o
	ranlib $(INCLDIR)/libcomparestructs_debug.a
	#
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/rational_mathematics.o $(SRCDIRLIB)/symlib/rational_mathematics.f90
	ar ru $(INCLDIR)/librational_debug.a $(DBGDIR_GT)/rational_mathematics.o $(DBGDIR_GT)/utilities.o
	ranlib $(INCLDIR)/librational_debug.a
	#
	$(FC) -c $(FFLAGS) -fpic $(DBGFFLAGS) -o $(DBGDIR_GT)/combinatorics.o $(SRCDIRLIB)/symlib/combinatorics.f90
	ar ru $(INCLDIR)/libcombinatorics_debug.a $(DBGDIR_GT)/combinatorics.o
	ranlib $(INCLDIR)/libcombinatorics_debug.a
	#
	$(FC) -c $(FFLAGS) -fpic  $(DBGFFLAGS) -I $(INCLDIR) -o $(DBGDIR_GT)/kpointgeneration.o $(SRCDIRLIB)/kgridlib/kpointgeneration.f90
	$(FC) -c $(FFLAGS) -fpic  $(DBGFFLAGS) -I $(INCLDIR) -o $(DBGDIR_GT)/wrap_kpoints.o $(SRCDIRLIB)/kgridlib/wrap_kpoints.f90
	ar ru $(INCLDIR)/libkptgen_debug.a $(DBGDIR_GT)/kpointgeneration.o $(DBGDIR_GT)/wrap_kpoints.o
	ranlib $(INCLDIR)/libkptgen_debug.a
	#
endif


#
# GTlib release rules-------------------------------------------------------------
RELGTlib:
ifeq ($(GT),T)
	#
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/num_types.o $(SRCDIRLIB)/symlib/num_types.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/utilities.o $(SRCDIRLIB)/symlib/utilities.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/numerical_utilities.o $(SRCDIRLIB)/symlib/numerical_utilities.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/vector_matrix_utilities.o $(SRCDIRLIB)/symlib/vector_matrix_utilities.f90
	ar ru $(INCLDIR)/libutils.a $(RELDIR_GT)/num_types.o $(RELDIR_GT)/utilities.o $(RELDIR_GT)/numerical_utilities.o $(RELDIR_GT)/vector_matrix_utilities.o
	ranlib $(INCLDIR)/libutils.a
	#
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/symmetry.o $(SRCDIRLIB)/symlib/symmetry.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/group_theory.o $(SRCDIRLIB)/symlib/group_theory.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/itertools.o $(SRCDIRLIB)/symlib/itertools.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/classes.o $(SRCDIRLIB)/symlib/classes.f90
	ar ru $(INCLDIR)/libsym.a $(RELDIR_GT)/symmetry.o $(RELDIR_GT)/group_theory.o $(RELDIR_GT)/itertools.o $(RELDIR_GT)/classes.o
	ranlib $(INCLDIR)/libsym.a
	#
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/compare_structures.o $(SRCDIRLIB)/symlib/compare_structures.f90
	ar ru $(INCLDIR)/libcomparestructs.a $(RELDIR_GT)/compare_structures.o
	ranlib $(INCLDIR)/libcomparestructs.a
	#
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/rational_mathematics.o $(SRCDIRLIB)/symlib/rational_mathematics.f90
	ar ru $(INCLDIR)/librational.a $(RELDIR_GT)/rational_mathematics.o $(RELDIR_GT)/utilities.o
	ranlib $(INCLDIR)/librational.a
	#
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -o $(RELDIR_GT)/combinatorics.o $(SRCDIRLIB)/symlib/combinatorics.f90
	ar ru $(INCLDIR)/libcombinatorics.a $(RELDIR_GT)/combinatorics.o
	ranlib $(INCLDIR)/libcombinatorics.a
	#
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -I $(INCLDIR) -o $(RELDIR_GT)/kpointgeneration.o $(SRCDIRLIB)/kgridlib/kpointgeneration.f90
	$(FC) -c $(FFLAGS) -fpic $(RELFFLAGS) -I $(INCLDIR) -o $(RELDIR_GT)/wrap_kpoints.o $(SRCDIRLIB)/kgridlib/wrap_kpoints.f90
	ar ru $(INCLDIR)/libkptgen.a $(RELDIR_GT)/kpointgeneration.o $(RELDIR_GT)/wrap_kpoints.o
	ranlib $(INCLDIR)/libkptgen.a
	#
endif


#
# test rules--------------------------------------------------------------------
test: prepDBG $(DBGEXE) moveTSTS


#
# Other rules-------------------------------------------------------------------
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
