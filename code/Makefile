CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -O3
CFLAGSDEBUG=-lm -lgsl -lgslcblas -Wall -g

test: test.o connect.o transport.o devices.o useful.o  mapping.o ../libs/matrices.o ../libs/hubbard.o ../libs/graphene_gf.o ../libs/suscept.o ../libs/static.o ../libs/greenfns.o ../libs/cubature.o
	$(CC)  test.o connect.o transport.o devices.o useful.o mapping.o ../libs/matrices.o ../libs/hubbard.o ../libs/graphene_gf.o ../libs/suscept.o ../libs/static.o ../libs/greenfns.o ../libs/cubature.o $(CFLAGS) -o test

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
	
../libs/matrices.o: ../libs/matrices.c ../libs/matrices.h
	$(CC) -c -g ../libs/matrices.c -o ../libs/matrices.o

../libs/hubbard.o: ../libs/hubbard.c ../libs/hubbard.h
	$(CC) -c ../libs/hubbard.c -o ../libs/hubbard.o

./libs/graphene_gf.o: ../libs/graphene_gf.c ../libs/graphene_gf.h
	$(CC) -c ../libs/graphene_gf.c -o ../libs/graphene_gf.o

../libs/suscept.o: ../libs/suscept.c ../libs/suscept.h 
	$(CC) -c ../libs/suscept.c -o ../libs/suscept.o

../libs/static.o: ../libs/static.c ../libs/static.h 
	$(CC) -c ../libs/static.c -o ../libs/static.o
	
../libs/greenfns.o: ../libs/greenfns.c ../libs/greenfns.h
	$(CC) -c ../libs/greenfns.c -o ../libs/greenfns.o
	
../libs/cubature.o: ../libs/cubature.c ../libs/cubature.h
	$(CC) -c ../libs/cubature.c -o ../libs/cubature.o

	
	
clean: FRC
	rm *.o ../libs/*.o
	
FRC: 		