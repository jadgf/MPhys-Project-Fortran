main = kmesh.f90
OBJmain = kmesh.o
FFLAGS = -O3 -ffree-line-length-0

 
FC = gfortran
LIB=-llapack -lblas

kmesh : ${OBJmain} 
	${FC} ${FFLAGS} -o kmesh kmesh.o ${LIB}

${OBJmain} : ${@:.o=.f90} ${main}
	${FC} ${FFLAGS} -c ${@:.o=.f90}

