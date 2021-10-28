CPP = icpc
CPPFLAG = -std=c++11 -qopenmp -mkl -DEIGEN_USE_MKL_ALL -O3 
EIGEN = ~/apps/Eigen3
GSL = ~/apps/gsl-2.6
TBLIS = ~/apps/tblis-1.2.0
MAIN = main.o mol.o scf.o ccsd.o molint.o intTrans.o
CHECK = check.o mol.o scf.o ccsd.o molint.o intTrans.o

test.exe: ${MAIN}
	${CPP} ${CPPFLAG} -I ${EIGEN} -I ${GSL} -I ${TBLIS}/include -L ${GSL}/.libs -l gsl -L ${TBLIS}/lib -l tblis ${MAIN} -o test.exe 

check.exe: ${CHECK}
	${CPP} ${CPPFLAG} -I ${EIGEN} -I ${GSL} -I ${TBLIS}/include -L ${GSL}/.libs -l gsl -L ${TBLIS}/lib -l tblis ${CHECK} -o check.exe 


%.o: %.cpp
	$(CPP) $(CPPFLAG) -c $< -o $@ -I ${EIGEN} -I ${GSL} -I ${TBLIS}/include  


clean:
	rm *.o *.exe
