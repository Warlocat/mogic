CPP = icpc
CPPFLAG = -O3 -std=c++11 -qopenmp -mkl -DEIGEN_USE_MKL_ALL
EIGEN = ~/apps/Eigen3
MAIN = main.o mol.o

test.exe: ${MAIN}
	${CPP} ${CPPFLAG} -I ${EIGEN} ${MAIN} -o test.exe


%.o: %.cpp
	$(CPP) $(CPPFLAG) -c $< -o $@ -I ${EIGEN}


clean:
	rm *.o *.exe
