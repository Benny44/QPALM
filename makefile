IDIR=include
ODIR=obj
SRCDIR=src
EDIR=examples
BDIR=lib

all: $(BDIR)/libqpalm.a demo

lib: $(BDIR)/libqpalm.a

CC=gcc
CFLAGS=-I$(IDIR) -ISuiteSparse/include -fPIC -O3 -DPROFILING 
LIBS=-lm 
CHOLMOD_LIBS=-lcholmod -lamd -lcolamd -lsuitesparseconfig
CHOLMOD_LIB_INCLUDES=-LSuiteSparse/lib

_DEPS = qpalm.h scaling.h util.h lin_alg.h validate.h linesearch.h types.h constants.h lbfgs.h global_opts.h termination.h cs.h cholmod_interface.h newton.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

__OBJ = qpalm.o scaling.o util.o lin_alg.o validate.o linesearch.o lbfgs.o termination.o cs.o cholmod_interface.o newton.o
_OBJ = $(patsubst %,$(ODIR)/%, $(__OBJ))
OBJ = $(EDIR)/qpalm_demo.o $(_OBJ)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(CHOLMOD_LIBS) $(CHOLMOD_LIB_INCLUDES)

$(EDIR)/%.o: $(EDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(CHOLMOD_LIBS) $(CHOLMOD_LIB_INCLUDES)

$(BDIR)/libqpalm.a: $(_OBJ)
	ar rcs $@ $^ 

demo: $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(CHOLMOD_LIBS) $(CHOLMOD_LIB_INCLUDES)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/*.a $(EDIR)/*.o demo
