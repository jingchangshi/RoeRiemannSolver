ccomp = g++
cflag = -Wall -O3 -std=c++11
main.out : main.cpp Solution.h
	$(ccomp) $(cflag) main.cpp -o main.out
run :
	./main.out
plot :
	python plot.py
clean :
	rm *.out *.txt *.dat *.png *.pyc
