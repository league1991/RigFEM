ifndef REDUCEDFORCEMODEL
REDUCEDFORCEMODEL=REDUCEDFORCEMODEL

ifndef CLEANFOLDER
CLEANFOLDER=REDUCEDFORCEMODEL
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
REDUCEDFORCEMODELOBJECTS=reducedForceModel.o reducedForceModelWithHessian.o

# the libraries this library depends on
REDUCEDFORCEMODELLIBS=

# the headers in this library
REDUCEDFORCEMODELHEADERS=reducedForceModel.h reducedForceModelWithHessian.h

REDUCEDFORCEMODELOBJECTS_FILENAMES=$(addprefix $(L)/reducedForceModel/, $(REDUCEDFORCEMODELOBJECTS))
REDUCEDFORCEMODELHEADER_FILENAMES=$(addprefix $(L)/reducedForceModel/, $(REDUCEDFORCEMODELHEADERS))
REDUCEDFORCEMODELLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(REDUCEDFORCEMODELLIBS))
REDUCEDFORCEMODELLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(REDUCEDFORCEMODELLIBS))

include $(REDUCEDFORCEMODELLIB_MAKEFILES)

all: $(L)/reducedForceModel/libreducedForceModel.a

$(L)/reducedForceModel/libreducedForceModel.a: $(REDUCEDFORCEMODELOBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/reducedForceModel/*.h $(L)/include

$(REDUCEDFORCEMODELOBJECTS_FILENAMES): %.o: %.cpp $(REDUCEDFORCEMODELLIB_FILENAMES) $(REDUCEDFORCEMODELHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), REDUCEDFORCEMODEL)
clean: cleanreducedForceModel
endif

deepclean: cleanreducedForceModel

cleanreducedForceModel:
	$(RM) $(REDUCEDFORCEMODELOBJECTS_FILENAMES) $(L)/reducedForceModel/libreducedForceModel.a

endif

