CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

.PHONY : all
all : rieszPyramid

Riesz_FFT: Riesz_FFT.o pgm.o fft.o
	gcc $(LDFLAGS) $^ -o $@ -lm
	
rieszPyramid: rieszPyramid.o rieszPyramidFilters.o pgm.o
	gcc $(LDFLAGS) $^ -o $@ -lm
	
#exo_2D_riesz: exo_2D_riesz.o pgm.o fft.o
#	gcc $(LDFLAGS) $^ -o $@ -lm

#exo_2D_gradient: exo_2D_gradient.o pgm.o fft.o
#	gcc $(LDFLAGS) $^ -o $@ -lm

%.o: %.c
	gcc $(CFLAGS) $< 

.PHONY : clean	
clean:
	rm -f *.o
	rm -f Riesz_FFT
	rm -f rieszPyramid
	rm -f *.pgm
	rm -f *~
#	rm -f exo_2D_riesz
#   rm -f exo_2D_gradient

