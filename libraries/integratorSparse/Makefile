ifndef INTEGRATORSPARSE
INTEGRATORSPARSE=INTEGRATORSPARSE

ifndef CLEANFOLDER
CLEANFOLDER=INTEGRATORSPARSE
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# NOTE: unlike the other library makefiles, this makefile adds $(PARDISO_INCLUDE) to the compiler invocation

# the object files to be compiled for this library
INTEGRATORSPARSEOBJECTS=centralDifferencesSparse.o eulerSparse.o implicitBackwardEulerSparse.o implicitNewmarkSparse.o integratorBaseSparse.o

# the libraries this library depends on
INTEGRATORSPARSELIBS=integrator performanceCounter insertRows sparseSolver forceModel

# the headers in this library
INTEGRATORSPARSEHEADERS=centralDifferencesSparse.h eulerSparse.h implicitBackwardEulerSparse.h implicitNewmarkSparse.h integratorBaseSparse.h

INTEGRATORSPARSEOBJECTS_FILENAMES=$(addprefix $(L)/integratorSparse/, $(INTEGRATORSPARSEOBJECTS))
INTEGRATORSPARSEHEADER_FILENAMES=$(addprefix $(L)/integratorSparse/, $(INTEGRATORSPARSEHEADERS))
INTEGRATORSPARSELIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(INTEGRATORSPARSELIBS))
INTEGRATORSPARSELIB_FILENAMES=$(call GET_LIB_FILENAMES, $(INTEGRATORSPARSELIBS))

include $(INTEGRATORSPARSELIB_MAKEFILES)

all: $(L)/integratorSparse/libintegratorSparse.a

$(L)/integratorSparse/libintegratorSparse.a: $(INTEGRATORSPARSEOBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/integratorSparse/*.h $(L)/include

$(INTEGRATORSPARSEOBJECTS_FILENAMES): %.o: %.cpp $(INTEGRATORSPARSELIB_FILENAMES) $(INTEGRATORSPARSEHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), INTEGRATORSPARSE)
clean: cleanintegratorSparse
endif

deepclean: cleanintegratorSparse

cleanintegratorSparse:
	$(RM) $(INTEGRATORSPARSEOBJECTS_FILENAMES) $(L)/integratorSparse/libintegratorSparse.a

endif

