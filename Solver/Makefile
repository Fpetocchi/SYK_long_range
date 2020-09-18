PALM_ROOT = $(HOME)/work/ALPS
include $(PALM_ROOT)/share/alps/include.mk

ECHO = echo
RM = rm -rf

.C.o:
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@  $<


all: sw

sw: impurity.o main.o
	$(CXX) -o $@ $(CPPFLAGS) $(CXXFLAGS) main.o impurity.o $(LDFLAGS) $(LIBS)

MPI_sw: impurity.o main.o
	$(CXX) -o $@ $(CPPFLAGS_MPI) $(CXXFLAGS) main.o impurity.o $(LDFLAGS_MPI)         $(LIBS_MPI)

clean:
	$(RM) ti_files *.o *~  sw MPI_sw

impurity.o: impurity.C impurity.h

evaluate: evaluate.o
	$(CXX) -o $@ $(CPPFLAGS) $(CXXFLAGS) evaluate.o $(LDFLAGS) $(LIBS)
