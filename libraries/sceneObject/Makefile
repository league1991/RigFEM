ifndef SCENEOBJECT
SCENEOBJECT=SCENEOBJECT

ifndef CLEANFOLDER
CLEANFOLDER=SCENEOBJECT
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
SCENEOBJECT_OBJECTS=sceneObject.o sceneObjectWithRestPosition.o sceneObject6DOF.o sceneObjectDeformable.o sceneObjectDeformable6DOF.o
ifneq ($(CG_LIBS),)
SCENEOBJECT_OBJECTS+=sceneObjectDeformableGPU.o sceneObjectDeformableGPU6DOF.o
endif

# the libraries this library depends on
SCENEOBJECT_LIBS=objMesh lighting minivector
ifneq ($(CG_LIBS),)
SCENEOBJECT_LIBS+=objMeshGPUDeformer
endif

# the headers in this library
SCENEOBJECT_HEADERS=sceneObject.h sceneObjectWithRestPosition.h sceneObject6DOF.h sceneObjectDeformable.h sceneObjectDeformable6DOF.h
ifneq ($(CG_LIBS),)
SCENEOBJECT_HEADERS+=sceneObjectDeformableGPU.h sceneObjectDeformableGPU6DOF.h
endif

SCENEOBJECT_OBJECTS_FILENAMES=$(addprefix $(L)/sceneObject/, $(SCENEOBJECT_OBJECTS))
SCENEOBJECT_HEADER_FILENAMES=$(addprefix $(L)/sceneObject/, $(SCENEOBJECT_HEADERS))
SCENEOBJECT_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(SCENEOBJECT_LIBS))
SCENEOBJECT_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(SCENEOBJECT_LIBS))

include $(SCENEOBJECT_LIB_MAKEFILES)

all: $(L)/sceneObject/libsceneObject.a

$(L)/sceneObject/libsceneObject.a: $(SCENEOBJECT_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/sceneObject/*.h $(L)/include

$(SCENEOBJECT_OBJECTS_FILENAMES): %.o: %.cpp $(SCENEOBJECT_LIB_FILENAMES) $(SCENEOBJECT_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(CG_INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), SCENEOBJECT)
clean: cleansceneObject
endif

deepclean: cleansceneObject

cleansceneObject:
	$(RM) $(SCENEOBJECT_OBJECTS_FILENAMES) $(L)/sceneObject/libsceneObject.a

endif

