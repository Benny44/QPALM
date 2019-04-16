IDIR=include
ODIR=obj
SRCDIR=src
EDIR=examples
BDIR=lib
TESTDIR=tests

all: directories lib demo

test: 
	(cd $(TESTDIR) && $(MAKE) )

lib: $(BDIR)/libqpalm.a

ifndef CC
	CC=gcc
endif

CFLAGS=-I$(IDIR) -Isuitesparse/include -fPIC -DPROFILING -Wall -Wextra -DDLONG -fopenmp -fexceptions -fprofile-arcs -ftest-coverage
CHOLMOD_LIBS=-lcholmod -lamd -lcolamd -lsuitesparseconfig -lcamd -lccolamd -lmetis -lm
CHOLMOD_LIB_INCLUDE+=-Lsuitesparse/lib -Isuitesparse/metis-5.1.0/include

LIBS+=$(CHOLMOD_LIBS)
LDLIBS +=$(CHOLMOD_LIB_INCLUDE)

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

_DEPS = qpalm.h scaling.h util.h lin_alg.h validate.h linesearch.h types.h constants.h global_opts.h termination.h cholmod_interface.h newton.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

__OBJ = qpalm.o scaling.o util.o lin_alg.o validate.o linesearch.o termination.o cholmod_interface.o newton.o
_OBJ = $(patsubst %,$(ODIR)/%, $(__OBJ))
OBJ = $(EDIR)/qpalm_demo.o $(_OBJ)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS) $(LDLIBS) 

$(EDIR)/%.o: $(EDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS) $(LDLIBS)

$(BDIR)/libqpalm.a: $(_OBJ)
	ar rcs $@ $^ 

demo: $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDLIBS)

.PHONY: clean directories

clean:
	rm -f $(ODIR)/*.o $(BDIR)/*.a $(EDIR)/*.o demo
	(cd $(TESTDIR) && $(MAKE) clean)

directories:
	mkdir -p obj lib

