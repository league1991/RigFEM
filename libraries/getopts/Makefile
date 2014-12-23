ifndef GETOPTS
GETOPTS=GETOPTS

ifndef CLEANFOLDER
CLEANFOLDER=GETOPTS
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
GETOPTS_OBJECTS=getopts.o

# the libraries this library depends on
GETOPTS_LIBS=

# the headers in this library
GETOPTS_HEADERS=getopts.h


GETOPTS_OBJECTS_FILENAMES=$(addprefix $(L)/getopts/, $(GETOPTS_OBJECTS))
GETOPTS_HEADER_FILENAMES=$(addprefix $(L)/getopts/, $(GETOPTS_HEADERS))
GETOPTS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(GETOPTS_LIBS))
GETOPTS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(GETOPTS_LIBS))

include $(GETOPTS_LIB_MAKEFILES)

all: $(L)/getopts/libgetopts.a

$(L)/getopts/libgetopts.a: $(GETOPTS_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/getopts/*.h $(L)/include

$(GETOPTS_OBJECTS_FILENAMES): %.o: %.cpp $(GETOPTS_LIB_FILENAMES) $(GETOPTS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), GETOPTS)
clean: cleangetopts
endif

deepclean: cleangetopts

cleangetopts:
	$(RM) $(GETOPTS_OBJECTS_FILENAMES) $(L)/getopts/libgetopts.a

endif
