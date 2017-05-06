/* Pyramides gaussienne et laplacienne */

#include "pgm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


/**
 * Monte d'un étage dans la pyramide Gaussienne
 * @param   img = image au niveau l
 * @return  image de niveau l+1
 */
double** reduce(double** img, double* h, int size,  int nl_old, int nc_old, int nl, int nc) {
    double** res = alloue_image_double(nl,nc);
    for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	    for (int ky=0; ky<size; ky++) {
    	        for (int kx=0; kx<size; kx++) {
    	            int y = 2*j+ky-2;
                    if (y<0) {
    	                y = nc_old+y;
    	            } else if (y>=nc_old) {
    	                y = y-nc_old;
                    }
                    int x = 2*i+kx-2;
                    if (x<0) {
    	                x = nl_old+x;
    	            } else if (x>=nl_old) {
    	                x = x-nl_old;
                    } 	      
    	            res[i][j] += img[x][y]*h[kx]*h[ky];
                }
    	    }
    	}
    }
    return res;
}

/**
 * Expansion x4 d'une image de la pyramide de Gauss (par interpolation)
 */
double** expand(double** img, double* h, int size, int nl_old, int nc_old, int nl, int nc) {
    double** res = alloue_image_double(nl,nc);
    for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	    for (int n=0; n<size; n++) {
    	        for (int m=0; m<size; m++) {
    	            // Seulement si les indices sont des entiers
    	            if (((j-m+2) % 2 == 0) && ((i-n+2) % 2 == 0)) {
        	            int y = (j-m+2)/2;
                        if (y<0) {
        	                y = nc_old+y;
        	            } else if (y>=nc_old) {
        	                y = y-nc_old;
                        }
                        int x = (i-n+2)/2;
                        if (x<0) {
        	                x = nl_old+x;
        	            } else if (x>=nl_old) {
        	                x = x-nl_old;
                        } 	      
        	            res[i][j] += img[x][y]*h[m]*h[n];
    	            }
                }
    	    }
    	    res[i][j] = 4*res[i][j];
    	}
    }
    return res;
}






int main (int ac, char **av) {

    // Noyau de lissage
	int size = 5;
    double* h = (double*) calloc(size,sizeof(double));
    h[0] = h[4] = 0.05;
    h[1] = h[3] = 0.25;
    h[2] = 0.4;
    
    // Paramètres de l'image en entrée
    int N, Ml, Mc, Nb;
    int nl, nc;
    
	unsigned char **im=NULL;
	double **im1;

	
	// Vérification du nombre de paramètres
    if (ac < 6) {
        printf("Usage : %s img_input N Ml Mc Nb\n", av[0]); exit(1);
    }
    
    // Lecture des paramètres
    im = lectureimagepgm(av[1],&nl,&nc);
    if (im==NULL) {
        puts("Lecture image impossible");
        exit(1);
    }
    N = atoi(av[2]);
    Ml = atoi(av[3]);
    Mc = atoi(av[4]);
    Nb = atoi(av[5]);
    if (Nb > N) {
        puts("On ne peut pas avoir plus d'étage que N");
        exit(1);
    }

    // Conversion de l'image de unsigned char** en double**
    im1 = imuchar2double(im,nl,nc);
    
    // Calcul des deux pyramides
    int* dimL = (int*) calloc(Nb+1,sizeof(int));
    int* dimC = (int*) calloc(Nb+1,sizeof(int));
    for (int k=0; k<Nb+1; k++) {
        dimL[k] = (double)Ml*(double)pow(2,N-k)+1;
        dimC[k] = (double)Mc*(double)pow(2,N-k)+1;
    }
    
    if ((nl != dimL[0]) || (nc != dimC[0])) {
        puts("Paramètres invalides");
        exit(1);
    }

    double*** pyramidGauss = (double***) calloc(Nb,sizeof(double**));
    double*** pyramidLap = (double***) calloc(Nb,sizeof(double**));
    double*** reduceGauss = (double***) calloc(Nb+1,sizeof(double**));
    reduceGauss[0] = im1;

    
    for (int k=0; k<Nb; k++) {
        reduceGauss[k+1] = reduce(reduceGauss[k], h, size, dimL[k], dimC[k], dimL[k+1], dimC[k+1]);
        double*** expands = (double***) calloc(k+2,sizeof(double**));
        expands[k+1] = reduceGauss[k+1];
        for (int l=k; l>=0; l--) {     
            expands[l] = expand(expands[l+1], h, size, dimL[l+1], dimC[l+1], dimL[l], dimC[l]);
        }
        pyramidGauss[k] = expands[0];
        // Libérationd de expands sauf expands[0]
        for (int l=1; l<k+1; l++) {
             libere_image_double(expands[l]);
        }
        free(expands);
    }

    
    for (int k=0; k<Nb; k++) {
        pyramidLap[k] = alloue_image_double(nl,nc);
        for (int i=0; i<nl; i++) {
    	    for (int j=0; j<nc; j++) {
    	        if (k==0) {
    	            pyramidLap[k][i][j] = im1[i][j]-pyramidGauss[0][i][j];
                } else {
                    pyramidLap[k][i][j] = pyramidGauss[k-1][i][j]-pyramidGauss[k][i][j];
                }
            }
        }
    }

    
    // Pyramide de Gauss
    for (int k=0; k<Nb; k++) {
        char* indice = (char*) calloc(1,sizeof(char));
        sprintf(indice, "%i", k+1);
        char filename[20];
        strcpy(filename, "Gauss_level");
        strcat(filename, indice);
        strcat(filename, ".pgm");
        ecritureimagepgm(filename, imdouble2uchar(pyramidGauss[k],nl,nc), nl, nc);
    }
    
    // Pyramide de Laplace
    for (int k=0; k<Nb; k++) {
        char* indice = (char*) calloc(1,sizeof(char));
        sprintf(indice, "%i", k);
        char filename[20];
        strcpy(filename, "Laplace_level");
        strcat(filename, indice);
        strcat(filename, ".pgm");
        ecritureimagepgm(filename, scale(pyramidLap[k],nl,nc), nl, nc);
    }

	// TEST
	double** test = alloue_image_double(nl,nc);
	for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	    test[i][j] = pyramidLap[0][i][j]+pyramidLap[1][i][j]+pyramidLap[2][i][j]+pyramidGauss[2][i][j];
        }
    }
    ecritureimagepgm("TEST.pgm", imdouble2uchar(test,nl,nc), nl, nc);
	
	
	
	
    return 0;
}  
