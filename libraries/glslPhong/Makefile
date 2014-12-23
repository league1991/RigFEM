ifndef GLSLPHONG
GLSLPHONG=GLSLPHONG

ifndef CLEANFOLDER
CLEANFOLDER=GLSLPHONG
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
GLSLPHONG_OBJECTS=glslPhong.o

# the libraries this library depends on
GLSLPHONG_LIBS=

# the headers in this library
GLSLPHONG_HEADERS=glslPhong.h


GLSLPHONG_OBJECTS_FILENAMES=$(addprefix $(L)/glslPhong/, $(GLSLPHONG_OBJECTS))
GLSLPHONG_HEADER_FILENAMES=$(addprefix $(L)/glslPhong/, $(GLSLPHONG_HEADERS))
GLSLPHONG_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(GLSLPHONG_LIBS))
GLSLPHONG_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(GLSLPHONG_LIBS))

include $(GLSLPHONG_LIB_MAKEFILES)

all: $(L)/glslPhong/libglslPhong.a

$(L)/glslPhong/libglslPhong.a: $(GLSLPHONG_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/glslPhong/*.h $(L)/include

$(GLSLPHONG_OBJECTS_FILENAMES): %.o: %.cpp $(GLSLPHONG_LIB_FILENAMES) $(GLSLPHONG_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), GLSLPHONG)
clean: cleanglslPhong
endif

deepclean: cleanglslPhong

cleanglslPhong:
	$(RM) $(GLSLPHONG_OBJECTS_FILENAMES) $(L)/glslPhong/libglslPhong.a

endif

