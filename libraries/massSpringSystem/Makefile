ifndef MASSSPRINGSYSTEM
MASSSPRINGSYSTEM=MASSSPRINGSYSTEM

ifndef CLEANFOLDER
CLEANFOLDER=MASSSPRINGSYSTEM
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
MASSSPRINGSYSTEM_OBJECTS=massSpringSystemFromObjMeshConfigFile.o massSpringSystemFromObjMesh.o massSpringSystemFromTetMeshConfigFile.o massSpringSystemFromTetMesh.o massSpringSystemMT.o massSpringSystem.o renderSprings.o massSpringSystemFromCubicMesh.o massSpringSystemFromCubicMeshConfigFile.o

# the libraries this library depends on
MASSSPRINGSYSTEM_LIBS=objMesh volumetricMesh configFile sparseMatrix

# the headers in this library
MASSSPRINGSYSTEM_HEADERS=massSpringSystemFromObjMeshConfigFile.h massSpringSystemFromObjMesh.h massSpringSystemFromTetMeshConfigFile.h massSpringSystemFromTetMesh.h massSpringSystem.h massSpringSystemMT.h renderSprings.h massSpringSystemFromCubicMesh.h massSpringSystemFromCubicMeshConfigFile.h

MASSSPRINGSYSTEM_OBJECTS_FILENAMES=$(addprefix $(L)/massSpringSystem/, $(MASSSPRINGSYSTEM_OBJECTS))
MASSSPRINGSYSTEM_HEADER_FILENAMES=$(addprefix $(L)/massSpringSystem/, $(MASSSPRINGSYSTEM_HEADERS))
MASSSPRINGSYSTEM_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(MASSSPRINGSYSTEM_LIBS))
MASSSPRINGSYSTEM_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(MASSSPRINGSYSTEM_LIBS))

include $(MASSSPRINGSYSTEM_LIB_MAKEFILES)

all: $(L)/massSpringSystem/libmassSpringSystem.a

$(L)/massSpringSystem/libmassSpringSystem.a: $(MASSSPRINGSYSTEM_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/massSpringSystem/*.h $(L)/include

$(MASSSPRINGSYSTEM_OBJECTS_FILENAMES): %.o: %.cpp $(MASSSPRINGSYSTEM_LIB_FILENAMES) $(MASSSPRINGSYSTEM_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), MASSSPRINGSYSTEM)
clean: cleanmassSpringSystem
endif

deepclean: cleanmassSpringSystem

cleanmassSpringSystem:
	$(RM) $(MASSSPRINGSYSTEM_OBJECTS_FILENAMES) $(L)/massSpringSystem/libmassSpringSystem.a

endif
