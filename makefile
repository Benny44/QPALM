IDIR=include
ODIR=obj
SRCDIR=src
EDIR=examples
BDIR=lib

all: $(BDIR)/libqpalm.a demo

CC=gcc-8
CFLAGS=-I$(IDIR) -fPIC -O3 -DPROFILING
LIBS=-lm

_DEPS = qpalm.h scaling.h util.h lin_alg.h validate.h linesearch.h types.h constants.h lbfgs.h global_opts.h termination.h cs.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

__OBJ = qpalm.o scaling.o util.o lin_alg.o validate.o linesearch.o lbfgs.o termination.o cs.o
_OBJ = $(patsubst %,$(ODIR)/%, $(__OBJ))
OBJ = $(EDIR)/qpalm_demo.o $(_OBJ)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EDIR)/%.o: $(EDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/libqpalm.a: $(_OBJ)
	ar rcs $@ $^ 

demo: $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/*.a $(EDIR)/*.o demo
