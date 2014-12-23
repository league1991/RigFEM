ifndef LOADLIST
LOADLIST=LOADLIST

ifndef CLEANFOLDER
CLEANFOLDER=LOADLIST
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
LOADLIST_OBJECTS=loadList.o

# the libraries this library depends on
LOADLIST_LIBS=

# the headers in this library
LOADLIST_HEADERS=loadList.h

LOADLIST_OBJECTS_FILENAMES=$(addprefix $(L)/loadList/, $(LOADLIST_OBJECTS))
LOADLIST_HEADER_FILENAMES=$(addprefix $(L)/loadList/, $(LOADLIST_HEADERS))
LOADLIST_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(LOADLIST_LIBS))
LOADLIST_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(LOADLIST_LIBS))

include $(LOADLIST_LIB_MAKEFILES)

all: $(L)/loadList/libloadList.a

$(L)/loadList/libloadList.a: $(LOADLIST_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/loadList/*.h $(L)/include

$(LOADLIST_OBJECTS_FILENAMES): %.o: %.cpp $(LOADLIST_LIB_FILENAMES) $(LOADLIST_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), LOADLIST)
clean: cleanloadList
endif

deepclean: cleanloadList

cleanloadList:
	$(RM) $(LOADLIST_OBJECTS_FILENAMES) $(L)/loadList/libloadList.a

endif
