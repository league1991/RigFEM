ifndef REDUCEDDYNAMICSOLVERRT
REDUCEDDYNAMICSOLVERRT=REDUCEDDYNAMICSOLVERRT

ifndef CLEANFOLDER
CLEANFOLDER=REDUCEDDYNAMICSOLVERRT
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this utility
REDUCEDDYNAMICSOLVERRT_OBJECTS=reducedDynamicSolver-rt.o initGraphics.o

# the libraries this utility depends on
# note: hashTable and rigidBody are not needed for this utility; they are included on this list so that they are compiled 
REDUCEDDYNAMICSOLVERRT_LIBS=objMeshGPUDeformer objMesh imageIO modalMatrix sceneObjectReduced sceneObject objMeshGPUDeformer objMesh imageIO integratorSparse integratorDense integrator reducedStvk stvk reducedElasticForceModel reducedForceModel forceModel matrix matrixIO lighting configFile volumetricMesh loadList getopts camera minivector openGLHelper hashTable rigidBodyDynamics

# the headers in this library
REDUCEDDYNAMICSOLVERRT_HEADERS=initGraphics.h

REDUCEDDYNAMICSOLVERRT_LINK=$(addprefix -l, $(REDUCEDDYNAMICSOLVERRT_LIBS)) $(SPOOLES_LIB) $(BLASLAPACK_LIB) $(GLUI_LIB) $(CG_LIBS) $(IMAGE_LIBS) $(STANDARD_LIBS)

REDUCEDDYNAMICSOLVERRT_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/reducedDynamicSolver-rt/, $(REDUCEDDYNAMICSOLVERRT_OBJECTS))
REDUCEDDYNAMICSOLVERRT_HEADER_FILENAMES=$(addprefix $(R)/utilities/reducedDynamicSolver-rt/, $(REDUCEDDYNAMICSOLVERRT_HEADERS))
REDUCEDDYNAMICSOLVERRT_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(REDUCEDDYNAMICSOLVERRT_LIBS))
REDUCEDDYNAMICSOLVERRT_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(REDUCEDDYNAMICSOLVERRT_LIBS))

include $(REDUCEDDYNAMICSOLVERRT_LIB_MAKEFILES) ../../libraries/glui/Makefile

all: $(R)/utilities/reducedDynamicSolver-rt/reducedDynamicSolver-rt

CURRENT_DIR = $(shell pwd)
$(R)/utilities/reducedDynamicSolver-rt/reducedDynamicSolver-rt: $(REDUCEDDYNAMICSOLVERRT_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(REDUCEDDYNAMICSOLVERRT_OBJECTS) $(REDUCEDDYNAMICSOLVERRT_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(REDUCEDDYNAMICSOLVERRT_OBJECTS_FILENAMES): %.o: %.cpp $(REDUCEDDYNAMICSOLVERRT_LIB_FILENAMES) $(REDUCEDDYNAMICSOLVERRT_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(CG_INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), REDUCEDDYNAMICSOLVERRT)
clean: cleanreducedDynamicSolver-rt
endif

deepclean: cleanreducedDynamicSolver-rt

cleanreducedDynamicSolver-rt:
	$(RM) $(REDUCEDDYNAMICSOLVERRT_OBJECTS_FILENAMES) $(R)/utilities/reducedDynamicSolver-rt/reducedDynamicSolver-rt

endif

