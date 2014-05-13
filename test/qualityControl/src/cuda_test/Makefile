INCLUDES=-I/home/nsarma/lapack/lapack-3.5.0/lapacke/include -I/opt/cuda-toolkit/5.5.22/include

LIB=/home/nsarma/lapack/lapack-3.5.0/liblapacke.a /home/nsarma/lapack/lapack-3.5.0/liblapack.a /home/nsarma/lapack/lapack-3.5.0/librefblas.a 
LIB+=/usr/lib64/libgfortran.so.3 
LIB+=/opt/cuda-toolkit/5.5.22/lib64/libcudart.so.5.5.22 /opt/cuda-toolkit/5.5.22/lib64/libcublas.so.5.5.22

CC = gcc

test: test.o matrix_functions.o test_cuda.o test_utils.o
	$(CC) -o test test.o matrix_functions.o test_cuda.o test_utils.o -lm $(LIB)

test.o: test.c test.h define_t.h
	$(CC) -O3 -c test.c -o test.o

matrix_functions.o: matrix_functions.c matrix_functions.h
	$(CC) $(INCLUDES) -O3 -c matrix_functions.c -o matrix_functions.o

test_cuda.o: test_cuda.cu test_cuda.h define_t.h
	nvcc  $(INCLUDES) -c test_cuda.cu -o test_cuda.o

test_utils.o: test_utils.cpp test_utils.h define_t.h
	$(CC) -O3 -c test_utils.cpp -o test_utils.o

clean: 
	rm -f *.o test
