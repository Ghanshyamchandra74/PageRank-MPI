CC=mpiicpc --std=c++11
CS=icpc --std=c++11

all: bin/MPI bin/Sequential 

bin/MPI: bin/MPI.o
	$(CC) -o bin/MPI bin/MPI.o

bin/MPI.o: src/MPI.cpp
	$(CC) -o bin/MPI.o -c src/MPI.cpp

bin/Sequential: bin/Sequential.o
	$(CS) -o bin/Sequential bin/Sequential.o

bin/Sequential.o: src/Sequential.cpp
	$(CS) -o bin/Sequential.o -c src/Sequential.cpp

clean:
	rm -f bin/MPI_Scan.o bin/MPI_Scan
	rm -f bin/Sequential.o bin/Sequential
