FC	= f77
FFLAGS	= -O3

# Among the 3 possibilities below, uncomment the appropriate
# case for your environment and comment the others.

# case 1: works when LAPACK and BLAS are installed.
#OBJLIBS =
#LIBS    = -llapack -lblas

# case 2: works when LAPACK is not installed but BLAS is.
#LIBS    = -lblas
#OBJLIBS = lapack.o

# case 3: works when neither LAPACK nor BLAS are installed.
OBJLIBS = blas.o lapack.o
LIBS	=

#############
#############

EXPOKIT   = expokit.o
#EXPOKIT   = dgpadm.o dspadm.o zgpadm.o zhpadm.o                    \
#            dgchbv.o dschbv.o zgchbv.o          dnchbv.o znchbv.o  \
#            dgexpv.o dsexpv.o zgexpv.o zhexpv.o                    \
#            dgphiv.o dsphiv.o zgphiv.o zhphiv.o 

MATAID    = mataid.o
#MATAID    = dgmatv.o zgmatv.o \
#            dgcnvr.o zgcnvr.o loadhb.o

OBJS      = clock.o $(EXPOKIT) $(MATAID) $(OBJLIBS)

sample_g: $(OBJS) sample_g.o
	$(FC) -o sample_g sample_g.o $(OBJS) $(LIBS)

sample_b: $(OBJS) sample_b.o
	$(FC) -o sample_b sample_b.o $(OBJS) $(LIBS)

sample_z: $(OBJS) sample_z.o
	$(FC) -o sample_z sample_z.o $(OBJS) $(LIBS)

sample_m: $(OBJS) sample_m.o
	$(FC) -o sample_m sample_m.o $(OBJS) $(LIBS)
 
sample_d: $(OBJS) sample_d.o
	$(FC) -o sample_d sample_d.o $(OBJS) $(LIBS)

sample_p: $(OBJS) sample_p.o
	$(FC) -o sample_p sample_p.o $(OBJS) $(LIBS)

.f.o:;  $(FC) $(FFLAGS) -c $<
