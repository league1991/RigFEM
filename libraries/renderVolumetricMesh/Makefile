ifndef RENDERVOLUMETRICMESH
RENDERVOLUMETRICMESH=RENDERVOLUMETRICMESH

ifndef CLEANFOLDER
CLEANFOLDER=RENDERVOLUMETRICMESH
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
RENDERVOLUMETRICMESH_OBJECTS=renderVolumetricMesh.o

# the libraries this library depends on
RENDERVOLUMETRICMESH_LIBS=volumetricMesh openGLHelper

# the headers in this library
RENDERVOLUMETRICMESH_HEADERS=renderVolumetricMesh.h

RENDERVOLUMETRICMESH_OBJECTS_FILENAMES=$(addprefix $(L)/renderVolumetricMesh/, $(RENDERVOLUMETRICMESH_OBJECTS))
RENDERVOLUMETRICMESH_HEADER_FILENAMES=$(addprefix $(L)/renderVolumetricMesh/, $(RENDERVOLUMETRICMESH_HEADERS))
RENDERVOLUMETRICMESH_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(RENDERVOLUMETRICMESH_LIBS))
RENDERVOLUMETRICMESH_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(RENDERVOLUMETRICMESH_LIBS))

include $(RENDERVOLUMETRICMESH_LIB_MAKEFILES)

all: $(L)/renderVolumetricMesh/librenderVolumetricMesh.a

$(L)/renderVolumetricMesh/librenderVolumetricMesh.a: $(RENDERVOLUMETRICMESH_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/renderVolumetricMesh/*.h $(L)/include

$(RENDERVOLUMETRICMESH_OBJECTS_FILENAMES): %.o: %.cpp $(RENDERVOLUMETRICMESH_LIB_FILENAMES) $(RENDERVOLUMETRICMESH_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), RENDERVOLUMETRICMESH)
clean: cleanrenderVolumetricMesh
endif

deepclean: cleanrenderVolumetricMesh

cleanrenderVolumetricMesh:
	$(RM) $(RENDERVOLUMETRICMESH_OBJECTS_FILENAMES) $(L)/renderVolumetricMesh/librenderVolumetricMesh.a

endif
