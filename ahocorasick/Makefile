LIBNAME := libahocorasick.a
CFLAGS := -Wall
HEADERS := actypes.h ahocorasick.h node.h

$(LIBNAME): ahocorasick.o node.o
	ar -cvq $(LIBNAME) ahocorasick.o node.o

ahocorasick.o: ahocorasick.c $(HEADERS)
	g++ -c ahocorasick.c $(CFLAGS)

node.o: node.c $(HEADERS)
	g++ -c node.c $(CFLAGS)

clean:
	rm -f ahocorasick.o node.o $(LIBNAME)
