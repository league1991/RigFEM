ifndef INTEGRATORDENSE
INTEGRATORDENSE=INTEGRATORDENSE

ifndef CLEANFOLDER
CLEANFOLDER=INTEGRATORDENSE
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# NOTE: unlike the other library makefiles, this makefile adds $(PARDISO_INCLUDE) to the compiler invocation

# the object files to be compiled for this library
INTEGRATORDENSE_OBJECTS=implicitNewmarkDenseMulti1D.o integratorMulti1D.o centralDifferencesDense.o implicitBackwardEulerDense.o implicitNewmarkDense.o integratorBaseDense.o

# the libraries this library depends on
INTEGRATORDENSE_LIBS=matrix integrator performanceCounter reducedForceModel

# the headers in this library
INTEGRATORDENSE_HEADERS=implicitNewmarkDenseMulti1D.h integratorMulti1D.h centralDifferencesDense.h implicitBackwardEulerDense.h implicitNewmarkDense.h integratorBaseDense.h

INTEGRATORDENSE_OBJECTS_FILENAMES=$(addprefix $(L)/integratorDense/, $(INTEGRATORDENSE_OBJECTS))
INTEGRATORDENSE_HEADER_FILENAMES=$(addprefix $(L)/integratorDense/, $(INTEGRATORDENSE_HEADERS))
INTEGRATORDENSE_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(INTEGRATORDENSE_LIBS))
INTEGRATORDENSE_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(INTEGRATORDENSE_LIBS))

include $(INTEGRATORDENSE_LIB_MAKEFILES)

all: $(L)/integratorDense/libintegratorDense.a

$(L)/integratorDense/libintegratorDense.a: $(INTEGRATORDENSE_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/integratorDense/*.h $(L)/include

$(INTEGRATORDENSE_OBJECTS_FILENAMES): %.o: %.cpp $(INTEGRATORDENSE_LIB_FILENAMES) $(INTEGRATORDENSE_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), INTEGRATORDENSE)
clean: cleanintegratorDense
endif

deepclean: cleanintegratorDense

cleanintegratorDense:
	$(RM) $(INTEGRATORDENSE_OBJECTS_FILENAMES) $(L)/integratorDense/libintegratorDense.a

endif

