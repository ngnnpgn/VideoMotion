#ifndef _FFT
#define _FFT 1
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"


int ispowerof2(int x);



int fft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int dimx, int dimy);



int ifft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int dimx, int dimy);




void fftshift( double** imsr, double** imsi, double** imdr, double** imdi, int nl, int nc );




double** padimdforfft(double** im, int* pnl, int* pnc);



double** padimucforfft(unsigned char** im, int* pnl, int* pnc);
