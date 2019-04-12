IDIR=include
ODIR=obj
SRCDIR=src
EDIR=examples
BDIR=lib

all: directories lib tests

__OBJ = test_lin_alg.o
_OBJ = $(patsubst %,$(ODIR)/%, $(__OBJ))
OBJ = run_all_tests.o $(_OBJ)

_DEPS = test_lin_alg.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

ifndef CC
	CC=gcc
endif

CFLAGS=-I$(IDIR) -I../$(IDIR) -I../suitesparse/include -fPIC -O3 -Wall -Wextra -fopenmp -fexceptions `pkg-config --cflags cunit` `pkg-config --libs cunit` -fuse-ld=gold
CHOLMOD_LIBS=-lcholmod -lamd -lcolamd -lsuitesparseconfig -lcamd -lccolamd -lmetis -lm
CHOLMOD_LIB_INCLUDE+=-L../suitesparse/lib -I../suitesparse/metis-5.1.0/include

LIBS+=$(CHOLMOD_LIBS)
LDLIBS+=$(CHOLMOD_LIB_INCLUDE)


#We need blas and lapack to compile. The user can specify this by running make BLAS="-lblas_library -llapack_library" BLAS_PATH=path/to/blas
ifndef BLAS
	BLAS=-lmwblas -lmwlapack
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
	./run_all_tests
	

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o run_all_tests


