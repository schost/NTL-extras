
INC_DIR = ../include
LIB_DIR = ../lib

OBJS = $(wildcard $(LIB_DIR)/*)

SRC = $(wildcard test*cpp)
TESTS = $(patsubst %.cpp, u%, $(SRC))
EXEC = $(patsubst %.cpp, %, $(SRC))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lpthread -lgmp -lm

BENCH_MACHINE = eschost-desktop

.PHONY: clean

clean:
	rm -f $(EXEC)

all:	dat process plot

dat:	block mulvec mulmat threshold crossover

utest_%:	$(OBJS)
	if test -e test_$*.cpp; then $(CXX) -o test_$* test_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi

block:	utest_cauchy_like_block
	./test_cauchy_like_block 1 > block.FFT0.$(shell hostname).dat
	./test_cauchy_like_block 0 > block.65537.$(shell hostname).dat

mulvec:	utest_cauchy_like_mul_vector
	./test_cauchy_like_mul_vector 1 > mul_vector.FFT0.$(shell hostname).dat
	./test_cauchy_like_mul_vector 0 > mul_vector.65537.$(shell hostname).dat

mulmat:	utest_cauchy_like_mul_matrix
	./test_cauchy_like_mul_matrix 1 > mul_matrix.FFT0.$(shell hostname).dat
	./test_cauchy_like_mul_matrix 0 > mul_matrix.65537.$(shell hostname).dat

threshold:	utest_cauchy_like_find_threshold
	./test_cauchy_like_find_threshold 1 > threshold-dac.FFT0.$(shell hostname).dat
	./test_cauchy_like_find_threshold 0 > threshold-dac.65537.$(shell hostname).dat

crossover:	utest_cauchy_like_invert
	./test_cauchy_like_invert 1 > crossover.FFT0.$(shell hostname).dat
	./test_cauchy_like_invert 0 > crossover.65537.$(shell hostname).dat

process: 
	awk -f threshold-block.awk block.FFT0.$(BENCH_MACHINE).dat > threshold-block.FFT0.$(BENCH_MACHINE).dat
	awk -f threshold-block.awk block.65537.$(BENCH_MACHINE).dat > threshold-block.65537.$(BENCH_MACHINE).dat
	awk -f ratio-block.awk block.FFT0.$(BENCH_MACHINE).dat > ratio-block.FFT0.$(BENCH_MACHINE).dat
	awk -f ratio-block.awk block.65537.$(BENCH_MACHINE).dat > ratio-block.65537.$(BENCH_MACHINE).dat
	awk -f threshold-mul.awk mul_matrix.FFT0.$(BENCH_MACHINE).dat > threshold-mul-matrix.FFT0.$(BENCH_MACHINE).dat
	awk -f threshold-mul.awk mul_matrix.65537.$(BENCH_MACHINE).dat > threshold-mul-matrix.65537.$(BENCH_MACHINE).dat
	awk -f ratio-mul.awk mul_matrix.FFT0.$(BENCH_MACHINE).dat > ratio-mul-matrix.FFT0.$(BENCH_MACHINE).dat
	awk -f ratio-mul.awk mul_matrix.65537.$(BENCH_MACHINE).dat > ratio-mul-matrix.65537.$(BENCH_MACHINE).dat
	awk -f find_mininum.awk threshold-dac.65537.$(BENCH_MACHINE).dat > threshold-dac-alpha.65537.$(BENCH_MACHINE).dat
	awk -f find_mininum.awk threshold-dac.FFT0.$(BENCH_MACHINE).dat > threshold-dac-alpha.FFT0.$(BENCH_MACHINE).dat
	awk -f crossover-invert.awk crossover.FFT0.$(BENCH_MACHINE).dat > crossover-invert.FFT0.$(BENCH_MACHINE).dat
	awk -f crossover-invert.awk crossover.65537.$(BENCH_MACHINE).dat > crossover-invert.65537.$(BENCH_MACHINE).dat
	awk -f crossover-mul-vector.awk mul_vector.FFT0.$(BENCH_MACHINE).dat > crossover-mul-vector.FFT0.$(BENCH_MACHINE).dat
	awk -f crossover-mul-vector.awk mul_vector.65537.$(BENCH_MACHINE).dat > crossover-mul-vector.65537.$(BENCH_MACHINE).dat

plot:
	gnuplot -e "infile1='threshold-block.65537.$(BENCH_MACHINE).dat'" -e "infile2='threshold-block.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='threshold-block.$(BENCH_MACHINE).pdf'" threshold-block.plt
	gnuplot -e "infile1='ratio-block.65537.$(BENCH_MACHINE).dat'" -e "infile2='ratio-block.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='ratio-block.$(BENCH_MACHINE).pdf'" ratio-block.plt
	gnuplot -e "infile1='threshold-mul-matrix.65537.$(BENCH_MACHINE).dat'" -e "infile2='threshold-mul-matrix.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='threshold-mul.$(BENCH_MACHINE).pdf'" threshold-mul.plt
	gnuplot -e "infile1='ratio-mul-matrix.65537.$(BENCH_MACHINE).dat'" -e "infile2='ratio-mul-matrix.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='ratio-mul-matrix.$(BENCH_MACHINE).pdf'" ratio-mul-matrix.plt
	gnuplot -e "infile1='threshold-dac-alpha.65537.$(BENCH_MACHINE).dat'" -e "infile2='threshold-dac-alpha.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='threshold-dac-alpha.$(BENCH_MACHINE).pdf'" threshold-fast-slow.plt
	gnuplot -e "infile1='crossover-invert.65537.$(BENCH_MACHINE).dat'" -e "infile2='crossover-invert.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='crossover-invert.$(BENCH_MACHINE).pdf'" crossover-invert.plt
	gnuplot -e "infile1='crossover-mul-vector.65537.$(BENCH_MACHINE).dat'" -e "infile2='crossover-mul-vector.FFT0.$(BENCH_MACHINE).dat'" -e "outfile='crossover-mul-vector.$(BENCH_MACHINE).pdf'" crossover-mul-vector.plt
