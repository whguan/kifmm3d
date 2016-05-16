all: libfmm3d_mpi

include ../makefile.mpi

LIB_SRC = 	kernel3d_mpi.cpp knlmat3d_mpi.cpp dense3d_mpi.cpp let3d_mpi.cpp fmm3d_mpi.cpp matmgnt3d_mpi.cpp fmm3d_setup_mpi.cpp fmm3d_eval_mpi.cpp fmm3d_check_mpi.cpp 

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	tt0.cpp	tt1.cpp

DEP     = 	$(LIB_SRC:.cpp=.d)

libfmm3d_mpi: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libfmm3d_mpi.a $(LIB_OBJ)
	$(RANLIB) libfmm3d_mpi.a

tt0: libfmm3d_mpi tt0.o
	${CXX} -o tt0 tt0.o libfmm3d_mpi.a ../common/libcommon.a ${LDFLAGS}

tt1: libfmm3d_mpi tt1.o
	${CXX} -o tt1 tt1.o libfmm3d_mpi.a ../common/libcommon.a ${LDFLAGS}

petsc_test: petsc_test.o
	${CXX} -o petsc_test petsc_test.o libfmm3d_mpi.a ../common/libcommon.a ${LDFLAGS}

-include $(DEP)
