all: libfmm3d

include ../makefile.opt

LIB_SRC = 	kernel3d.cpp knlmat3d.cpp dense3d.cpp let3d.cpp fmm3d.cpp matmgnt3d.cpp fmm3d_setup.cpp fmm3d_eval.cpp fmm3d_check.cpp 

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	tt.cpp

MEX_SRC = 	fmm3d_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d)

libfmm3d: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libfmm3d.a $(LIB_OBJ)
	$(RANLIB) libfmm3d.a

tt: libfmm3d tt.o
	${CXX} -o tt tt.o libfmm3d.a ../common/libcommon.a ${LDFLAGS}

matlab: libfmm3d.a ${MEX_SRC}
	${MEX} fmm3d_mex.cpp ${INCLUDES} libfmm3d.a ../common/libcommon.a ${LIBS}
	mv *_mex.mex* ../matlab

-include $(DEP)

#	${MEX} fmm3d_mex.cpp ${INCLUDES} libfmm3d.a ../common/libcommon.a librfftw.a libfftw.a liblapack.a libblas.a



