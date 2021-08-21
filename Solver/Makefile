#
# Compiler flags
#
CC = gcc
CXX = mpicxx

CFLAGS = -std=c++11 -lstdc++ -lm #-std=c++0x -stdlib=libc++ -lgsl -lgslcblas
DBGCFLAGS = -O0 -Wall -Wextra -fmax-errors=1000 -pedantic -ggdb -D_verb
RELCFLAGS = -O3 -fpermissive -DNDEBUG

#
# Project files-----------------------------------------------------------------
INCLDIR = include
EXECDIR = bin
SRCDIR = src
SRCS = main.cpp
CFLAGS +=-I $(INCLDIR)/eigen
OBJS = $(SRCS:.cpp=.o)
EXE  = ct_hyb


#
# Debug build settings----------------------------------------------------------
DBGDIR = debug
DBGEXE = $(EXECDIR)/$(EXE)_debug
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
#DBGCFLAGS += -odir $(DBGDIR)


#
# Release build settings--------------------------------------------------------
RELDIR = release
RELEXE = $(EXECDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
#RELCFLAGS += -odir $(RELDIR)


#
# All rules---------------------------------------------------------------------
.PHONY: all clean debug prep release remake


#
# Default build-----------------------------------------------------------------
all: prepREL release


#
# Debug rules-------------------------------------------------------------------
debug: prepDBG $(DBGEXE)

$(DBGDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) $(DBGCFLAGS) -c $< -o $@

$(DBGEXE): $(DBGOBJS)
	$(CXX) $(CFLAGS) $(DBGCFLAGS) $^ -o $@



#
# Release rules-----------------------------------------------------------------
release: prepREL $(RELEXE)

$(RELDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) $(RELCFLAGS) -c $< -o $@

$(RELEXE): $(RELOBJS)
	$(CXX) $(CFLAGS) $(RELCFLAGS) $^ -o $@


#
# Other rules-------------------------------------------------------------------
prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(EXECDIR)

prepREL:
	@mkdir -p $(RELDIR) $(EXECDIR)

prepDBG:
	@mkdir -p $(DBGDIR) $(EXECDIR)

remake: clean all

clean:
	rm -rf $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS) $(DBGDIR) $(RELDIR) $(EXECDIR)

info:
	echo "tumadre"
