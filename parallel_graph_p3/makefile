MPI_LSTDFLG = -lstdc++ -llapack -lcblas -lm
MPI_INCLUDE = -I/usr/include/
MPI_LIB = -L/usr/lib/atlas/
MPI_OBJS = parallel_graph

all:	${MPI_OBJS}
	rm -f *.o

parallel_graph.o: parallel_graph.cpp
	mpic++ -g -c parallel_graph.cpp -o parallel_graph.o ${MPI_INCLUDE}

parallel_graph: parallel_graph.o
	mpic++ parallel_graph.o -o parallel_graph ${MPI_LIB} ${MPI_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MPI_OBJS}