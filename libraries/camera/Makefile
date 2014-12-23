ifndef CAMERA
CAMERA=CAMERA

ifndef CLEANFOLDER
CLEANFOLDER=CAMERA
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
CAMERA_OBJECTS=camera.o

# the libraries this library depends on
CAMERA_LIBS=

# the headers in this library
CAMERA_HEADERS=camera.h


CAMERA_OBJECTS_FILENAMES=$(addprefix $(L)/camera/, $(CAMERA_OBJECTS))
CAMERA_HEADER_FILENAMES=$(addprefix $(L)/camera/, $(CAMERA_HEADERS))
CAMERA_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(CAMERA_LIBS))
CAMERA_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(CAMERA_LIBS))

include $(CAMERA_LIB_MAKEFILES)

all: $(L)/camera/libcamera.a

$(L)/camera/libcamera.a: $(CAMERA_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/camera/*.h $(L)/include

$(CAMERA_OBJECTS_FILENAMES): %.o: %.cpp $(CAMERA_LIB_FILENAMES) $(CAMERA_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), CAMERA)
clean: cleancamera
endif

deepclean: cleancamera

cleancamera:
	$(RM) $(CAMERA_OBJECTS_FILENAMES) $(L)/camera/libcamera.a

endif
