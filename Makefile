CPP = icpc
CPPFLAG = -std=c++11 -qopenmp -mkl -DEIGEN_USE_MKL_ALL -O3 -ipo
EIGEN = ~/apps/Eigen3
GSL = ~/apps/gsl-2.6
MAIN = main.o mol.o scf.o ccsd.o molint.o intTrans.o

test.exe: ${MAIN}
	${CPP} ${CPPFLAG} -I ${EIGEN} -I ${GSL} -L ${GSL}/.libs -l gsl ${MAIN} -o test.exe 


%.o: %.cpp
	$(CPP) $(CPPFLAG) -c $< -o $@ -I ${EIGEN} 


clean:
	rm *.o *.exe
