ifndef OPENGLHELPER
OPENGLHELPER=OPENGLHELPER

ifndef CLEANFOLDER
CLEANFOLDER=OPENGLHELPER
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
OPENGLHELPER_OBJECTS=openGLHelper.o printBitmap.o saveScreenShot.o

# the libraries this library depends on
OPENGLHELPER_LIBS=matrixIO imageIO

# the headers in this library
OPENGLHELPER_HEADERS=openGLHelper.h printBitmap.h saveScreenShot.h

OPENGLHELPER_OBJECTS_FILENAMES=$(addprefix $(L)/openGLHelper/, $(OPENGLHELPER_OBJECTS))
OPENGLHELPER_HEADER_FILENAMES=$(addprefix $(L)/openGLHelper/, $(OPENGLHELPER_HEADERS))
OPENGLHELPER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(OPENGLHELPER_LIBS))
OPENGLHELPER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(OPENGLHELPER_LIBS))

include $(OPENGLHELPER_LIB_MAKEFILES)

all: $(L)/openGLHelper/libopenGLHelper.a

$(L)/openGLHelper/libopenGLHelper.a: $(OPENGLHELPER_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/openGLHelper/*.h $(L)/include

$(OPENGLHELPER_OBJECTS_FILENAMES): %.o: %.cpp $(OPENGLHELPER_LIB_FILENAMES) $(OPENGLHELPER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), OPENGLHELPER)
clean: cleanopenGLHelper
endif

deepclean: cleanopenGLHelper

cleanopenGLHelper:
	$(RM) $(OPENGLHELPER_OBJECTS_FILENAMES) $(L)/openGLHelper/libopenGLHelper.a

endif
