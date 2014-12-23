ifndef MODALMATRIX
MODALMATRIX=MODALMATRIX

ifndef CLEANFOLDER
CLEANFOLDER=MODALMATRIX
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
MODALMATRIX_OBJECTS=modalMatrix.o

# the libraries this library depends on
MODALMATRIX_LIBS=

# the headers in this library
MODALMATRIX_HEADERS=modalMatrix.h

MODALMATRIX_OBJECTS_FILENAMES=$(addprefix $(L)/modalMatrix/, $(MODALMATRIX_OBJECTS))
MODALMATRIX_HEADER_FILENAMES=$(addprefix $(L)/modalMatrix/, $(MODALMATRIX_HEADERS))
MODALMATRIX_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(MODALMATRIX_LIBS))
MODALMATRIX_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(MODALMATRIX_LIBS))

include $(MODALMATRIX_LIB_MAKEFILES)

all: $(L)/modalMatrix/libmodalMatrix.a

$(L)/modalMatrix/libmodalMatrix.a: $(MODALMATRIX_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/modalMatrix/*.h $(L)/include

$(MODALMATRIX_OBJECTS_FILENAMES): %.o: %.cpp $(MODALMATRIX_LIB_FILENAMES) $(MODALMATRIX_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), MODALMATRIX)
clean: cleanmodalMatrix
endif

deepclean: cleanmodalMatrix

cleanmodalMatrix:
	$(RM) $(MODALMATRIX_OBJECTS_FILENAMES) $(L)/modalMatrix/libmodalMatrix.a

endif

