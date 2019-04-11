IDIR=include
ODIR=obj
SRCDIR=src
EDIR=examples
BDIR=lib

all: lib demo

lib: $(BDIR)/libqpalm.a

CC=gcc
CFLAGS=-I$(IDIR) -Isuitesparse/include -fPIC -O3 -DPROFILING -Wall -Wextra -DDLONG -fopenmp -fexceptions
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

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/*.a $(EDIR)/*.o demo
