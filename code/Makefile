CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -O3
CFLAGSDEBUG=-lm -lgsl -lgslcblas -Wall -g
MYLIBDIR=../libs

test: test.o connect.o useful_hops.o transport.o devices.o disorder.o analytic.o useful.o mapping.o matrices.o hubbard.o graphene_gf.o suscept.o static.o greenfns.o cubature.o 
	$(CC)  test.o connect.o useful_hops.o transport.o devices.o disorder.o analytic.o useful.o mapping.o matrices.o hubbard.o graphene_gf.o suscept.o static.o greenfns.o cubature.o  $(CFLAGS) -o test

test.o: test.c test.h
	$(CC) -c -g -O3 test.c

devices.o: devices.c devices.h
	$(CC) -c -g -O3 devices.c
	
useful.o: useful.c useful.h
	$(CC) -c -g -O3 useful.c
	
connect.o: connect.c connect.h
	$(CC) -c -g -O3 connect.c
	
transport.o: transport.c transport.h
	$(CC) -c -g -O3 transport.c

mapping.o: mapping.c mapping.h
	$(CC) -c -g -O3 mapping.c
	
useful_hops.o: useful_hops.c useful_hops.h
	$(CC) -c -g -O3 useful_hops.c
	
disorder.o: disorder.c disorder.h
	$(CC) -c -g -O3 disorder.c
	
analytic.o: analytic.c analytic.h 
	$(CC) -c -g -O3 analytic.c
	
matrices.o: $(MYLIBDIR)/matrices.c $(MYLIBDIR)/matrices.h
	$(CC) -c -g $(MYLIBDIR)/matrices.c -o matrices.o

hubbard.o: $(MYLIBDIR)/hubbard.c $(MYLIBDIR)/hubbard.h
	$(CC) -c $(MYLIBDIR)/hubbard.c -o hubbard.o

graphene_gf.o: $(MYLIBDIR)/graphene_gf.c $(MYLIBDIR)/graphene_gf.h
	$(CC) -c $(MYLIBDIR)/graphene_gf.c -o graphene_gf.o

suscept.o: $(MYLIBDIR)/suscept.c $(MYLIBDIR)/suscept.h 
	$(CC) -c $(MYLIBDIR)/suscept.c -o suscept.o

static.o: $(MYLIBDIR)/static.c $(MYLIBDIR)/static.h 
	$(CC) -c $(MYLIBDIR)/static.c -o static.o
	
greenfns.o: $(MYLIBDIR)/greenfns.c $(MYLIBDIR)/greenfns.h
	$(CC) -c $(MYLIBDIR)/greenfns.c -o greenfns.o
	
cubature.o: $(MYLIBDIR)/cubature.c $(MYLIBDIR)/cubature.h
	$(CC) -c $(MYLIBDIR)/cubature.c -o cubature.o

	
	
clean: FRC
	
FRC: 		