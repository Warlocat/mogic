CPP = icpc
CPPFLAG = -std=c++11 -qopenmp -mkl -DEIGEN_USE_MKL_ALL -O3 
EIGEN = ~/apps/Eigen3
MAIN = main.o mol.o scf.o ccsd.o

test.exe: ${MAIN}
	${CPP} ${CPPFLAG} -I ${EIGEN} ${MAIN} -o test.exe


%.o: %.cpp
	$(CPP) $(CPPFLAG) -c $< -o $@ -I ${EIGEN}


clean:
	rm *.o *.exe
