ifndef COROTATIONALLINEARFEM
COROTATIONALLINEARFEM=COROTATIONALLINEARFEM

ifndef CLEANFOLDER
CLEANFOLDER=COROTATIONALLINEARFEM
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
COROTATIONALLINEARFEM_OBJECTS=corotationalLinearFEM.o corotationalLinearFEMMT.o

# the libraries this library depends on
COROTATIONALLINEARFEM_LIBS=polarDecomposition volumetricMesh sparseMatrix

# the headers in this library
COROTATIONALLINEARFEM_HEADERS=corotationalLinearFEM.h corotationalLinearFEMMT.h


COROTATIONALLINEARFEM_OBJECTS_FILENAMES=$(addprefix $(L)/corotationalLinearFEM/, $(COROTATIONALLINEARFEM_OBJECTS))
COROTATIONALLINEARFEM_HEADER_FILENAMES=$(addprefix $(L)/corotationalLinearFEM/, $(COROTATIONALLINEARFEM_HEADERS))
COROTATIONALLINEARFEM_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(COROTATIONALLINEARFEM_LIBS))
COROTATIONALLINEARFEM_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(COROTATIONALLINEARFEM_LIBS))

include $(COROTATIONALLINEARFEM_LIB_MAKEFILES)

all: $(L)/corotationalLinearFEM/libcorotationalLinearFEM.a

$(L)/corotationalLinearFEM/libcorotationalLinearFEM.a: $(COROTATIONALLINEARFEM_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/corotationalLinearFEM/*.h $(L)/include

$(COROTATIONALLINEARFEM_OBJECTS_FILENAMES): %.o: %.cpp $(COROTATIONALLINEARFEM_LIB_FILENAMES) $(COROTATIONALLINEARFEM_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), COROTATIONALLINEARFEM)
clean: cleancorotationalLinearFEM
endif

deepclean: cleancorotationalLinearFEM

cleancorotationalLinearFEM:
	$(RM) $(COROTATIONALLINEARFEM_OBJECTS_FILENAMES) $(L)/corotationalLinearFEM/libcorotationalLinearFEM.a

endif
