CXX	= g++
CXXOBJS = main.o opt.o simplex.o simulated_annealing.o genetic_algorithms.o tree.o interior_point_method.o chaos_optimization.o optimize.o minimum_spanning_tree.o branch_and_cut.o

INCL	= Makefile
CXXINCL = opt.h simplex.h simulated_annealing.h genetic_algorithms.h	position.h	params.h	tree.h	utils.h interior_point_method.h function.h chaos_optimization.h minimum_spanning_tree.h branch_and_cut.h
CXXFLAGS = -O3 -Wall 

DEST	= /usr/local/bin
LIBS	= -L/usr/lib 
LIBS	+= -lm  

EXE	= Optimize 


EXE: ${CXXOBJS}
	${CXX} ${CXXFLAGS} ${CXXOBJS} ${LIBS} -o ${EXE} 

CXXOBJS: ${INCL} ${CXXINCL}

clean:	
	rm -f *.o *.gnu ${CXXOBJS} ${EXE}
