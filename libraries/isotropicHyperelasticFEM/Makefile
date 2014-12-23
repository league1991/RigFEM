ifndef IHFEM
IHFEM=IHFEM

ifndef CLEANFOLDER
CLEANFOLDER=IHFEM
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this library
IHFEM_OBJECTS=isotropicMaterial.o isotropicMaterialWithCompressionResistance.o MooneyRivlinIsotropicMaterial.o neoHookeanIsotropicMaterial.o StVKIsotropicMaterial.o homogeneousMooneyRivlinIsotropicMaterial.o homogeneousStVKIsotropicMaterial.o homogeneousNeoHookeanIsotropicMaterial.o isotropicHyperelasticFEM.o isotropicHyperelasticFEMMT.o

# the libraries this library depends on
IHFEM_LIBS=minivector volumetricMesh sparseMatrix

# the headers in this library
IHFEM_HEADERS=isotropicMaterial.h isotropicMaterialWithCompressionResistance.h MooneyRivlinIsotropicMaterial.h neoHookeanIsotropicMaterial.h StVKIsotropicMaterial.h homogeneousMooneyRivlinIsotropicMaterial.h homogeneousStVKIsotropicMaterial.h homogeneousNeoHookeanIsotropicMaterial.h isotropicHyperelasticFEM.h isotropicHyperelasticFEMMT.h

IHFEM_OBJECTS_FILENAMES=$(addprefix $(L)/isotropicHyperelasticFEM/, $(IHFEM_OBJECTS))
IHFEM_HEADER_FILENAMES=$(addprefix $(L)/isotropicHyperelasticFEM/, $(IHFEM_HEADERS))
IHFEM_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(IHFEM_LIBS))
IHFEM_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(IHFEM_LIBS))

include $(IHFEM_LIB_MAKEFILES)

all: $(L)/isotropicHyperelasticFEM/libisotropicHyperelasticFEM.a

$(L)/isotropicHyperelasticFEM/libisotropicHyperelasticFEM.a: $(IHFEM_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/isotropicHyperelasticFEM/*.h $(L)/include

$(IHFEM_OBJECTS_FILENAMES): %.o: %.cpp $(IHFEM_LIB_FILENAMES) $(IHFEM_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), IHFEM)
clean: cleanisotropicHyperelasticFEM
endif

deepclean: cleanisotropicHyperelasticFEM

cleanisotropicHyperelasticFEM:
	$(RM) $(IHFEM_OBJECTS_FILENAMES) $(L)/isotropicHyperelasticFEM/libisotropicHyperelasticFEM.a

endif

