all: libcommon

include ../makefile.mpi

LIB_SRC = 	vecmatop.cpp	svdrep.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

DEP     = 	$(LIB_SRC:.cpp=.d)

libcommon: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libcommon.a $(LIB_OBJ)
	$(RANLIB) libcommon.a

-include $(DEP)

