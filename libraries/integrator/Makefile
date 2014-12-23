ifndef INTEGRATOR
INTEGRATOR=INTEGRATOR

ifndef CLEANFOLDER
CLEANFOLDER=INTEGRATOR
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# NOTE: unlike the other library makefiles, this makefile adds $(PARDISO_INCLUDE) to the compiler invocation

# the object files to be compiled for this library
INTEGRATOROBJECTS=integratorBase.o getIntegratorSolver.o

# the libraries this library depends on
INTEGRATORLIBS=

# the headers in this library
INTEGRATORHEADERS=integratorBase.h getIntegratorSolver.h integratorSolverSelection.h

INTEGRATOROBJECTS_FILENAMES=$(addprefix $(L)/integrator/, $(INTEGRATOROBJECTS))
INTEGRATORHEADER_FILENAMES=$(addprefix $(L)/integrator/, $(INTEGRATORHEADERS))
INTEGRATORLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(INTEGRATORLIBS))
INTEGRATORLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(INTEGRATORLIBS))

include $(INTEGRATORLIB_MAKEFILES)

all: $(L)/integrator/libintegrator.a

$(L)/integrator/libintegrator.a: $(INTEGRATOROBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/integrator/*.h $(L)/include

$(INTEGRATOROBJECTS_FILENAMES): %.o: %.cpp $(INTEGRATORLIB_FILENAMES) $(INTEGRATORHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), INTEGRATOR)
clean: cleanintegrator
endif

deepclean: cleanintegrator

cleanintegrator:
	$(RM) $(INTEGRATOROBJECTS_FILENAMES) $(L)/integrator/libintegrator.a

endif

