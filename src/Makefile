FC = gfortran
FCFLAGS = -ggdb -c \
	  -fdefault-real-8 \
	  -Wall -Wsurprising -Wextra -Wunderflow -pedantic \
	  -fcheck=all \
	  -fbacktrace \
	  -std=f2003

FCFLAGS_FALSE = -c -g -fdefault-real-8 

FLFLAGS = -ggdb -fbacktrace

main: main.o fortran_gnuplot.o
	$(FC) $(FLFLAGS) -o main *.o

main.o: main.f03 fortran_gnuplot.mod repo.mod
	$(FC) $(FCFLAGS) main.f03

fortran_gnuplot.o: fortran_gnuplot.f03
	$(FC) $(FCFLAGS_FALSE) fortran_gnuplot.f03 

fortran_gnuplot.mod: fortran_gnuplot.f03
	$(FC) $(FCFLAGS_FALSE) fortran_gnuplot.f03

repo.o: repo.f03
	$(FC) $(FCFLAGS) repo.f03

repo.mod: repo.f03
	$(FC) $(FCFLAGS) repo.f03

run: main
	./main

clean:
	rm *.o *.mod *.dat

cleanall:
	rm main *.o *.mod 


