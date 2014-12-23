ifndef REDUCEDSTVK
REDUCEDSTVK=REDUCEDSTVK

ifndef CLEANFOLDER
CLEANFOLDER=REDUCEDSTVK
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
REDUCEDSTVKOBJECTS=StVKReducedInternalForces.o StVKReducedInternalForcesMT.o StVKReducedStiffnessMatrix.o StVKReducedHessianTensor.o StVKReducedInternalForcesMultiLoad.o

# the libraries this library depends on
REDUCEDSTVKLIBS=stvk modalMatrix matrix

# the headers in this library
REDUCEDSTVKHEADERS=StVKReducedInternalForces.h StVKReducedInternalForcesMT.h StVKReducedStiffnessMatrix.h StVKReducedHessianTensor.h StVKReducedInternalForcesMultiLoad.h

REDUCEDSTVKOBJECTS_FILENAMES=$(addprefix $(L)/reducedStvk/, $(REDUCEDSTVKOBJECTS))
REDUCEDSTVKHEADER_FILENAMES=$(addprefix $(L)/reducedStvk/, $(REDUCEDSTVKHEADERS))
REDUCEDSTVKLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(REDUCEDSTVKLIBS))
REDUCEDSTVKLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(REDUCEDSTVKLIBS))

include $(REDUCEDSTVKLIB_MAKEFILES)

all: $(L)/reducedStvk/libreducedStvk.a

$(L)/reducedStvk/libreducedStvk.a: $(REDUCEDSTVKOBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/reducedStvk/*.h $(L)/include

$(REDUCEDSTVKOBJECTS_FILENAMES): %.o: %.cpp $(REDUCEDSTVKLIB_FILENAMES) $(REDUCEDSTVKHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), REDUCEDSTVK)
clean: cleanreducedStvk
endif

deepclean: cleanreducedStvk

cleanreducedStvk:
	$(RM) $(REDUCEDSTVKOBJECTS_FILENAMES) $(L)/reducedStvk/libreducedStvk.a

endif

