ifndef REDUCEDELASTICFORCEMODEL
REDUCEDELASTICFORCEMODEL=REDUCEDELASTICFORCEMODEL

ifndef CLEANFOLDER
CLEANFOLDER=REDUCEDELASTICFORCEMODEL
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
REDUCEDELASTICFORCEMODELOBJECTS=reducedLinearForceModel.o reducedLinearStVKForceModel.o reducedMassSpringSystemForceModel.o reducedMassSpringSystemForceModelWithHessian.o reducedSpringForceModel.o reducedStVKForceModel.o reducedStVKForceModelWithHessian.o reducedSubspaceStVKForceModel.o

# the libraries this library depends on
REDUCEDELASTICFORCEMODELLIBS=reducedForceModel reducedStvk stvk massSpringSystem

# the headers in this library
REDUCEDELASTICFORCEMODELHEADERS=reducedLinearForceModel.h reducedLinearStVKForceModel.h reducedMassSpringSystemForceModel.h reducedMassSpringSystemForceModelWithHessian.h reducedSpringForceModel.h reducedStVKForceModel.h reducedStVKForceModelWithHessian.h reducedSubspaceStVKForceModel.h

REDUCEDELASTICFORCEMODELOBJECTS_FILENAMES=$(addprefix $(L)/reducedElasticForceModel/, $(REDUCEDELASTICFORCEMODELOBJECTS))
REDUCEDELASTICFORCEMODELHEADER_FILENAMES=$(addprefix $(L)/reducedElasticForceModel/, $(REDUCEDELASTICFORCEMODELHEADERS))
REDUCEDELASTICFORCEMODELLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(REDUCEDELASTICFORCEMODELLIBS))
REDUCEDELASTICFORCEMODELLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(REDUCEDELASTICFORCEMODELLIBS))

include $(REDUCEDELASTICFORCEMODELLIB_MAKEFILES)

all: $(L)/reducedElasticForceModel/libreducedElasticForceModel.a

$(L)/reducedElasticForceModel/libreducedElasticForceModel.a: $(REDUCEDELASTICFORCEMODELOBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/reducedElasticForceModel/*.h $(L)/include

$(REDUCEDELASTICFORCEMODELOBJECTS_FILENAMES): %.o: %.cpp $(REDUCEDELASTICFORCEMODELLIB_FILENAMES) $(REDUCEDELASTICFORCEMODELHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), REDUCEDELASTICFORCEMODEL)
clean: cleanreducedElasticForceModel
endif

deepclean: cleanreducedElasticForceModel

cleanreducedElasticForceModel:
	$(RM) $(REDUCEDELASTICFORCEMODELOBJECTS_FILENAMES) $(L)/reducedElasticForceModel/libreducedElasticForceModel.a

endif

