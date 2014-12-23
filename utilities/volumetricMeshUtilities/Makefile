ifndef VOLUTILS
VOLUTILS=VOLUTILS

ifndef CLEANFOLDER
CLEANFOLDER=VOLUTILS
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this utility
VOLUTILS_OBJECTS=

# the libraries this utility depends on
VOLUTILS_LIBS=volumetricMesh objMesh imageIO matrixIO getopts loadList sparseMatrix graph minivector

# the headers in this utility
VOLUTILS_HEADERS=

VOLUTILS_LINK=$(addprefix -l, $(VOLUTILS_LIBS)) $(SPOOLES_LIB) $(PARDISO_LIB) $(STANDARD_LIBS)


VOLUTILS_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/volumetricMeshUtilities/, $(VOLUTILS_OBJECTS))
VOLUTILS_HEADER_FILENAMES=$(addprefix $(R)/utilities/volumetricMeshUtilities/, $(VOLUTILS_HEADERS))
VOLUTILS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(VOLUTILS_LIBS))
VOLUTILS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(VOLUTILS_LIBS))

include $(VOLUTILS_LIB_MAKEFILES)

all: $(R)/utilities/volumetricMeshUtilities/generateMassMatrix $(R)/utilities/volumetricMeshUtilities/generateInterpolant $(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix $(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh

$(R)/utilities/volumetricMeshUtilities/generateMassMatrix: $(R)/utilities/volumetricMeshUtilities/generateMassMatrix.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(VOLUTILS_OBJECTS) $(R)/utilities/volumetricMeshUtilities/generateMassMatrix.cpp $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/volumetricMeshUtilities/generateInterpolant: $(R)/utilities/volumetricMeshUtilities/generateInterpolant.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(VOLUTILS_OBJECTS) $(R)/utilities/volumetricMeshUtilities/generateInterpolant.cpp $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix: $(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(VOLUTILS_OBJECTS) $(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix.cpp $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh: $(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(VOLUTILS_OBJECTS) $(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh.cpp $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(VOLUTILS_OBJECTS_FILENAMES): %.o: %.cpp $(VOLUTILS_LIB_FILENAMES) $(VOLUTILS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), VOLUTILS)
clean: cleanvolumetricMeshUtilities
endif

deepclean: cleanvolumetricMeshUtilities

cleanvolumetricMeshUtilities:
	$(RM) $(VOLUTILS_OBJECTS_FILENAMES) $(R)/utilities/volumetricMeshUtilities/generateMassMatrix $(R)/utilities/volumetricMeshUtilities/generateInterpolant $(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix $(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh

endif
