ifndef INSERTROWS
INSERTROWS=INSERTROWS

ifndef CLEANFOLDER
CLEANFOLDER=INSERTROWS
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
INSERTROWS_OBJECTS=insertRows.o

# the libraries this library depends on
INSERTROWS_LIBS=

# the headers in this library
INSERTROWS_HEADERS=insertRows.h


INSERTROWS_OBJECTS_FILENAMES=$(addprefix $(L)/insertRows/, $(INSERTROWS_OBJECTS))
INSERTROWS_HEADER_FILENAMES=$(addprefix $(L)/insertRows/, $(INSERTROWS_HEADERS))
INSERTROWS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(INSERTROWS_LIBS))
INSERTROWS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(INSERTROWS_LIBS))

include $(INSERTROWS_LIB_MAKEFILES)

all: $(L)/insertRows/libinsertRows.a

$(L)/insertRows/libinsertRows.a: $(INSERTROWS_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/insertRows/*.h $(L)/include

$(INSERTROWS_OBJECTS_FILENAMES): %.o: %.cpp $(INSERTROWS_LIB_FILENAMES) $(INSERTROWS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), INSERTROWS)
clean: cleaninsertRows
endif

deepclean: cleaninsertRows

cleaninsertRows:
	$(RM) $(INSERTROWS_OBJECTS_FILENAMES) $(L)/insertRows/libinsertRows.a

endif
