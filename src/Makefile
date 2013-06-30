IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -std=c99 -lm -O3
BIN = ../bin
ODIR=obj
LDIR =../lib

_DEPS = mmio.h functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = mmio.o functions.o preconditioning.o gmres.o cg.o main.o  
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BIN)/main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
