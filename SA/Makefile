CC = g++
IDIR = INCLUDE
ODIR = OBJS
CFLAGS = -O3 -Wall -I$(IDIR) -g -std=c++14

_DEPS = SA.h utils.h node.h model.h 

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = Main.o  SAopt.o
             
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

EXE = SA

$(ODIR)/%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: 
	all clean

all:
	$(MAKE) $(EXE)

$(EXE): $(OBJ) $(DEPS)
	$(CC) -o $(EXE) $(OBJ) $(CFLAGS) -lm

clean:
	/bin/rm -f $(ODIR)/*.o $(EXE) *~ ._* $(IDIR)/*~ $(IDIR)/._* 

