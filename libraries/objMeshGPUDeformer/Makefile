ifndef OBJMESHGPUDEFORMER
OBJMESHGPUDEFORMER=OBJMESHGPUDEFORMER

ifndef CLEANFOLDER
CLEANFOLDER=OBJMESHGPUDEFORMER
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
OBJMESHGPUDEFORMER_OBJECTS=objMeshGPUDeformer.o objMeshGPUDeformer_coarseToFine.o objMeshGPUDeformer_coarseToFine_fbo.o objMeshGPUDeformer_uUq.o objMeshGPUDeformer_uUq_fbo.o objMeshGPUDeformer_uUq_pbuffer.o objMeshGPUDeformer_uUq_setLighting.o vbo.o

# the libraries this library depends on
OBJMESHGPUDEFORMER_LIBS=objMesh lighting

# the headers in this library
OBJMESHGPUDEFORMER_HEADERS=objMeshGPUDeformer.h objMeshGPUDeformer_coarseToFine.h objMeshGPUDeformer_coarseToFine_fbo.h objMeshGPUDeformer_uUq.h objMeshGPUDeformer_uUq_fbo.h objMeshGPUDeformer_uUq_pbuffer.h objMeshGPUDeformer_uUq_setLighting.h vbo.h

OBJMESHGPUDEFORMER_OBJECTS_FILENAMES=$(addprefix $(L)/objMeshGPUDeformer/, $(OBJMESHGPUDEFORMER_OBJECTS))
OBJMESHGPUDEFORMER_HEADER_FILENAMES=$(addprefix $(L)/objMeshGPUDeformer/, $(OBJMESHGPUDEFORMER_HEADERS))
OBJMESHGPUDEFORMER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(OBJMESHGPUDEFORMER_LIBS))
OBJMESHGPUDEFORMER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(OBJMESHGPUDEFORMER_LIBS))

include $(OBJMESHGPUDEFORMER_LIB_MAKEFILES)

all: $(L)/objMeshGPUDeformer/libobjMeshGPUDeformer.a

$(L)/objMeshGPUDeformer/libobjMeshGPUDeformer.a: $(OBJMESHGPUDEFORMER_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/objMeshGPUDeformer/*.h $(L)/include

$(OBJMESHGPUDEFORMER_OBJECTS_FILENAMES): %.o: %.cpp $(OBJMESHGPUDEFORMER_LIB_FILENAMES) $(OBJMESHGPUDEFORMER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(CG_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), OBJMESHGPUDEFORMER)
clean: cleanobjMeshGPUDeformer
endif

deepclean: cleanobjMeshGPUDeformer

cleanobjMeshGPUDeformer:
	$(RM) $(OBJMESHGPUDEFORMER_OBJECTS_FILENAMES) $(L)/objMeshGPUDeformer/libobjMeshGPUDeformer.a

endif
