ifndef LIGHTING
LIGHTING=LIGHTING

ifndef CLEANFOLDER
CLEANFOLDER=LIGHTING
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
LIGHTING_OBJECTS=lighting.o

# the libraries this library depends on
LIGHTING_LIBS=configFile

# the headers in this library
LIGHTING_HEADERS=lighting.h

LIGHTING_OBJECTS_FILENAMES=$(addprefix $(L)/lighting/, $(LIGHTING_OBJECTS))
LIGHTING_HEADER_FILENAMES=$(addprefix $(L)/lighting/, $(LIGHTING_HEADERS))
LIGHTING_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(LIGHTING_LIBS))
LIGHTING_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(LIGHTING_LIBS))

include $(LIGHTING_LIB_MAKEFILES)

all: $(L)/lighting/liblighting.a

$(L)/lighting/liblighting.a: $(LIGHTING_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/lighting/*.h $(L)/include

$(LIGHTING_OBJECTS_FILENAMES): %.o: %.cpp $(LIGHTING_LIB_FILENAMES) $(LIGHTING_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), LIGHTING)
clean: cleanlighting
endif

deepclean: cleanlighting

cleanlighting:
	$(RM) $(LIGHTING_OBJECTS_FILENAMES) $(L)/lighting/liblighting.a

endif
