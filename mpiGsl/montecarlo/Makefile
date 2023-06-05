
CFLAGS = `pkg-config --cflags gsl` -O3
LFLAGS = `pkg-config --libs gsl`

all: seq para

para: para.c
	mpicc $(CFLAGS) para.c $(LFLAGS) -o para

seq: seq.c
	gcc $(CFLAGS) seq.c $(LFLAGS) -o seq

clean:
	rm -f para seq daten.dat *~

