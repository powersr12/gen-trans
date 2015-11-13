
CC=icc
DEF_MOD=/zdata/groups/common/nicpa/$(NPA_SOFTWAREVERSION)/$(NPA_NODETYPE)
GSL_PATH=$(DEF_MOD)/gsl/1.16/intel-15.0.1/lib
CFLAGS= -L$(GSL_PATH) -Wl,-rpath=$(GSL_PATH) -lgsl -mkl=parallel -lmkl_blas95_lp64 -O3
CFLAGSDEBUG=-lm -lgsl -lgslcblas -Wall -g

MYLIBDIR=../libs

test: test.o connect.o transport.o devices.o useful.o  mapping.o matrices.o hubbard.o graphene_gf.o suscept.o static.o greenfns.o cubature.o
	$(CC)  test.o connect.o transport.o devices.o useful.o mapping.o matrices.o hubbard.o graphene_gf.o suscept.o static.o greenfns.o cubature.o $(CFLAGS) -o test

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
	rm *.o
	
FRC: 		