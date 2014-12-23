ifndef DISPLAYOBJ
DISPLAYOBJ=DISPLAYOBJ

ifndef CLEANFOLDER
CLEANFOLDER=DISPLAYOBJ
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for these utilities
DISPLAYOBJ_OBJECTS=

# the libraries these utilities depend on
DISPLAYOBJ_LIBS=getopts objMesh imageIO glslPhong camera matrixIO lighting configFile loadList

# the headers in this library
DISPLAYOBJ_HEADERS=

DISPLAYOBJ_LINK=$(addprefix -l, $(DISPLAYOBJ_LIBS)) $(SPOOLES_LIB) $(PARDISO_LIB) $(GLEW_LIB) $(STANDARD_LIBS) $(GLUI_LIB)

DISPLAYOBJ_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/displayObj/, $(DISPLAYOBJ_OBJECTS))
DISPLAYOBJ_HEADER_FILENAMES=$(addprefix $(R)/utilities/displayObj/, $(DISPLAYOBJ_HEADERS))
DISPLAYOBJ_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(DISPLAYOBJ_LIBS))
DISPLAYOBJ_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(DISPLAYOBJ_LIBS))

include $(DISPLAYOBJ_LIB_MAKEFILES)

all: $(R)/utilities/displayObj/displayObj $(R)/utilities/displayObj/objMergeFiles

$(R)/utilities/displayObj/displayObj: $(R)/utilities/displayObj/displayObj.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(GLUI_INCLUDE) $(DISPLAYOBJ_OBJECTS) $^ $(DISPLAYOBJ_LINK) -Wl,-rpath,$(GLUI_DIR)/lib $(IMAGE_LIBS) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/displayObj/objMergeFiles: $(R)/utilities/displayObj/objMergeFiles.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(GLUI_INCLUDE) $(DISPLAYOBJ_OBJECTS) $^ $(DISPLAYOBJ_LINK) -Wl,-rpath,$(GLUI_DIR)/lib $(IMAGE_LIBS) -o $@; cp $@ $(R)/utilities/bin/

$(DISPLAYOBJ_OBJECTS_FILENAMES): %.o: %.cpp $(DISPLAYOBJ_LIB_FILENAMES) $(DISPLAYOBJ_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(IMAGE_INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), DISPLAYOBJ)
clean: cleandisplayObj
endif

deepclean: cleandisplayObj

cleandisplayObj:
	$(RM) $(R)/utilities/displayObj/displayObj 

endif

