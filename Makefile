all: seq par

seq:
	mpicc -Wall -o seq sequential.c

par:
	mpicc -Wall -o par parallel.c

