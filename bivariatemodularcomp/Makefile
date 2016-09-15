
test: test.o bivariatemodularcomp.o mat_ZZX.o mat_ZZ_pX.o
	g++ -std=c++14 -g test.o mat_ZZX.o bivariatemodularcomp.o mat_ZZ_pX.o -o test -lntl -lgmp -lm

test.o: test.cc
	g++ -std=c++14 -g test.cc -c

bivariatemodularcomp.o: bivariatemodularcomp.h bivariatemodularcomp.cc
	g++ -std=c++14 -g bivariatemodularcomp.cc -c

mat_ZZX.o: mat_ZZX.h mat_ZZX.cc
	g++ -std=c++14 -g mat_ZZX.cc -c

mat_ZZ_pX.o: mat_ZZ_pX.h mat_ZZ_pX.cc
	g++ -std=c++14 -g mat_ZZ_pX.cc -c

.PHONY: clean
clean:
	rm test *.o
