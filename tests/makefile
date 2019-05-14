IDIR=include
ODIR=obj
SRCDIR=src
EDIR=examples
BDIR=lib

all: directories lib tests

__OBJ = test_lin_alg.o test_cholmod_interface.o test_basic_qp.o test_prim_inf_qp.o test_dua_inf_qp.o test_degen_hess.o test_update.o test_validate.o test_nonconvex_qp.o
_OBJ = $(patsubst %,$(ODIR)/%, $(__OBJ))
OBJ = run_all_tests.o $(_OBJ)

_DEPS = test_lin_alg.h test_cholmod_interface.h test_basic_qp.h test_prim_inf_qp.h test_dua_inf_qp.h test_degen_hess.h test_update.h test_validate.h test_nonconvex_qp.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

ifndef CC
	CC=gcc
endif

CFLAGS=-I$(IDIR) -I../$(IDIR) -I../suitesparse/include -fPIC -Wall -DDLONG -Wextra -fopenmp -fexceptions 
CHOLMOD_LIBS=-lcholmod -lamd -lcolamd -lsuitesparseconfig -lcamd -lccolamd -lmetis -lm
CHOLMOD_LIB_INCLUDE+=-L../suitesparse/lib -I../suitesparse/metis-5.1.0/include

LIBS+=$(CHOLMOD_LIBS)
LDLIBS+=$(CHOLMOD_LIB_INCLUDE)

#Testing and checking coverage (also used in travis)
ifeq ($(COV),yes)
	CFLAGS+= `pkg-config --cflags cunit` `pkg-config --libs cunit` -fuse-ld=gold -fprofile-arcs -ftest-coverage
	BLAS=-lblas -llapack
endif

#We need blas and lapack to compile. The user can specify this by running make BLAS="-lblas_library -llapack_library" BLAS_PATH=path/to/blas
ifndef BLAS
	BLAS=-lblas -llapack
endif
LIBS+=$(BLAS)

ifdef BLAS_PATH
	BLAS_INCLUDE=-L$(BLAS_PATH)
	LDLIBS+=$(BLAS_INCLUDE)
endif


ifdef LAPACK_PATH
	LAPACK_INCLUDE=-L$(LAPACK_PATH)
	LDLIBS+=$(LAPACK_INCLUDE)
endif




directories:
	mkdir -p obj

lib: 
	(cd .. && $(MAKE) lib)

LIBS+=-lqpalm
LDLIBS+=-L../lib

tests: build_tests run_all_tests

build_tests: $(OBJ)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS) $(LDLIBS)

run_all_tests.o: run_all_tests.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS) $(LDLIBS)

run_all_tests: $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDLIBS)
	rm -f run_all_tests.o
	#For travis LD_LIBRARY_PATH seems to be required
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):../suitesparse/lib:../lib ./run_all_tests
	

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o run_all_tests


