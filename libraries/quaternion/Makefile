ifndef QUATERNION
QUATERNION=QUATERNION

ifndef CLEANFOLDER
CLEANFOLDER=QUATERNION
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
QUATERNIONOBJECTS=quaternion.o

# the libraries this library depends on
QUATERNIONLIBS=

# the headers in this library
QUATERNIONHEADERS=quaternion.h

QUATERNIONOBJECTS_FILENAMES=$(addprefix $(L)/quaternion/, $(QUATERNIONOBJECTS))
QUATERNIONHEADER_FILENAMES=$(addprefix $(L)/quaternion/, $(QUATERNIONHEADERS))
QUATERNIONLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(QUATERNIONLIBS))
QUATERNIONLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(QUATERNIONLIBS))

include $(QUATERNIONLIB_MAKEFILES)

all: $(L)/quaternion/libquaternion.a

$(L)/quaternion/libquaternion.a: $(QUATERNIONOBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/quaternion/*.h $(L)/include

$(QUATERNIONOBJECTS_FILENAMES): %.o: %.cpp $(QUATERNIONLIB_FILENAMES) $(QUATERNIONHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(BLASLAPACK_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), QUATERNION)
clean: cleanquaternion
endif

deepclean: cleanquaternion

cleanquaternion:
	$(RM) $(QUATERNIONOBJECTS_FILENAMES) $(L)/quaternion/libquaternion.a

endif

