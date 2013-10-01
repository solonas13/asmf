MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -msse3 -O3 -fomit-frame-pointer -funroll-loops -I./ahocorasick 
 
LFLAGS= -I $(HOME)/sdsl/include -L $(HOME)/sdsl/lib -lsdsl -ldivsufsort -ldivsufsort64 -lm -lahocorasick -L./ahocorasick -lz -ldatrie
 
EXE=    example
 
SRC=    example.cc sm.cc csm.cc lcp.cc aca.cc filter.cc
 
HD=     lfs.h sm.h csm.h lcp.h aca.h filter.h Makefile
 
# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cc .o 
 
OBJ=    $(SRC:.cc=.o) 
 
.cc.o: 
	$(CC) $(CFLAGS)-c $(LFLAGS) $< 
 
all:    $(EXE) 
 
$(EXE): $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS) 
 
$(OBJ): $(MF) $(HD) 
 
clean: 
	rm -f $(OBJ) $(EXE) *~
