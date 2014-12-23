ifndef POLARDECOMPOSITION
POLARDECOMPOSITION=POLARDECOMPOSITION

ifndef CLEANFOLDER
CLEANFOLDER=POLARDECOMPOSITION
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
POLARDECOMPOSITION_OBJECTS=polarDecomposition.o polarDecompositionGradient.o

# the libraries this library depends on
POLARDECOMPOSITION_LIBS=

# the headers in this library
POLARDECOMPOSITION_HEADERS=polarDecompositionGradient.h polarDecomposition.h

POLARDECOMPOSITION_OBJECTS_FILENAMES=$(addprefix $(L)/polarDecomposition/, $(POLARDECOMPOSITION_OBJECTS))
POLARDECOMPOSITION_HEADER_FILENAMES=$(addprefix $(L)/polarDecomposition/, $(POLARDECOMPOSITION_HEADERS))
POLARDECOMPOSITION_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(POLARDECOMPOSITION_LIBS))
POLARDECOMPOSITION_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(POLARDECOMPOSITION_LIBS))

include $(POLARDECOMPOSITION_LIB_MAKEFILES)

all: $(L)/polarDecomposition/libpolarDecomposition.a

$(L)/polarDecomposition/libpolarDecomposition.a: $(POLARDECOMPOSITION_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/polarDecomposition/*.h $(L)/include

$(POLARDECOMPOSITION_OBJECTS_FILENAMES): %.o: %.cpp $(POLARDECOMPOSITION_LIB_FILENAMES) $(POLARDECOMPOSITION_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), POLARDECOMPOSITION)
clean: cleanpolarDecomposition
endif

deepclean: cleanpolarDecomposition

cleanpolarDecomposition:
	$(RM) $(POLARDECOMPOSITION_OBJECTS_FILENAMES) $(L)/polarDecomposition/libpolarDecomposition.a

endif
