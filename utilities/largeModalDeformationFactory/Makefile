ifndef LARGEMODALDEFORMATIONFACTORY
LARGEMODALDEFORMATIONFACTORY=LARGEMODALDEFORMATIONFACTORY

ifndef CLEANFOLDER
CLEANFOLDER=LARGEMODALDEFORMATIONFACTORY
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this utility
LARGEMODALDEFORMATIONFACTORY_OBJECTS=cubicPolynomials.o fixedVertices.o frequencies.o largeModalDeformationFactory.o linearModes.o main.o nonlinearModes.o view.o renderingMesh.o simulationMesh.o canvas.o modalDerivatives.o sketch.o interpolate.o convert.o runtime.o StVKReducedInternalForcesWX.o

# the libraries this utility depends on
LARGEMODALDEFORMATIONFACTORY_LIBS=reducedStvk stvk reducedElasticForceModel reducedForceModel forceModel renderVolumetricMesh sparseSolver sparseMatrix volumetricMesh objMesh imageIO modalMatrix matrix matrixIO getopts insertRows loadList camera minivector openGLHelper

# the headers in this library
LARGEMODALDEFORMATIONFACTORY_HEADERS=StVKReducedInternalForcesWX.h canvas.h largeModalDeformationFactory.h states.h

CG_LIBS=-framework Cg
LARGEMODALDEFORMATIONFACTORY_LINK=$(addprefix -l, $(LARGEMODALDEFORMATIONFACTORY_LIBS)) $(ARPACK_LIB) $(SPOOLES_LIB) $(BLASLAPACK_LIB) $(PARDISO_LIB) $(FORTRAN_LIB) $(CG_LIBS) $(IMAGE_LIBS) $(STANDARD_LIBS)

LARGEMODALDEFORMATIONFACTORY_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/largeModalDeformationFactory/, $(LARGEMODALDEFORMATIONFACTORY_OBJECTS))
LARGEMODALDEFORMATIONFACTORY_HEADER_FILENAMES=$(addprefix $(R)/utilities/largeModalDeformationFactory/, $(LARGEMODALDEFORMATIONFACTORY_HEADERS))
LARGEMODALDEFORMATIONFACTORY_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(LARGEMODALDEFORMATIONFACTORY_LIBS))
LARGEMODALDEFORMATIONFACTORY_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(LARGEMODALDEFORMATIONFACTORY_LIBS))

include $(LARGEMODALDEFORMATIONFACTORY_LIB_MAKEFILES)

TARGET=largeModalDeformationFactory
all: $(TARGET) bundle

$(TARGET): $(LARGEMODALDEFORMATIONFACTORY_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(LARGEMODALDEFORMATIONFACTORY_OBJECTS) $(BIND_AT_LOAD) $(LARGEMODALDEFORMATIONFACTORY_LINK) `$(WX_CONFIG) --cxxflags --libs core,base,gl` -o $@; cp $@ $(R)/utilities/bin/

$(LARGEMODALDEFORMATIONFACTORY_OBJECTS_FILENAMES): %.o: %.cpp $(LARGEMODALDEFORMATIONFACTORY_LIB_FILENAMES) $(LARGEMODALDEFORMATIONFACTORY_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) `$(WX_CONFIG) --cxxflags` -c $(BLASLAPACK_INCLUDE) $(PARDISO_INCLUDE) $(INCLUDE) $< -o $@

bundle: $(TARGET) $(WX_DIR)/src/osx/carbon/Info.plist.in
	mkdir -p $(TARGET).app/Contents
	mkdir -p $(TARGET).app/Contents/MacOS
	mkdir -p $(TARGET).app/Contents/Resources
	
	sed -e "s/IDENTIFIER/`echo $(WX_DIR) | sed -e 's,\.\./,,g' | sed -e 's,/,.,g'`/" \
	-e "s/EXECUTABLE/$(TARGET)/" \
	-e "s/VERSION/$(WX_VERSION)/" \
	$(WX_DIR)/src/osx/carbon/Info.plist.in >$(TARGET).app/Contents/Info.plist
	
	echo -n "APPL????" >$(TARGET).app/Contents/PkgInfo
	
	ln -f $(TARGET)$(EXEEXT) $(TARGET).app/Contents/MacOS/$(TARGET)
	
	cp -f largeModalDeformationFactory.icns $(TARGET).app/Contents/Resources

ifeq ($(CLEANFOLDER), LARGEMODALDEFORMATIONFACTORY)
clean: cleanlargeModalDeformationFactory
endif

deepclean: cleanlargeModalDeformationFactory

cleanlargeModalDeformationFactory:
	$(RM) $(LARGEMODALDEFORMATIONFACTORY_OBJECTS_FILENAMES) $(R)/utilities/largeModalDeformationFactory/largeModalDeformationFactory

endif

