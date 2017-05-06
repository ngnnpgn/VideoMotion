

#include "pgm.h"
#include "rieszPyramidFilters.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * Application d'un filtre
 * @param   img		image à filtrer
 * @param	h		filtre de taille size*size
 * @param	nl,nc	dimensions
 * @return  image filtrée de même taille que img
 */
double** applyHighFilter(double** img, double** h, int nl, int nc) {
    double** res = alloue_image_double(nl,nc);
    for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	    for (int n=-4; n<=4; n++) {
    	        for (int m=-4; m<4; m++) {
    	        	int x = i+n;
    	        	int y = j+m;
                    if (x<0) {
    	                x = x+nl;
    	            } else if (x>=nl) {
    	                x = x-nl;
                    } 	   
                    if (y<0) {
    	                y = y+nc;
    	            } else if (y>=nc) {
    	                y = y-nc;
                    }
    	            res[i][j] += img[x][y]*h[n+4][m+4];
                }
    	    }
    	}
    }
    return res;
}


/**
 * Pyramide passe bas
 * @param   img = image au niveau l
 * @return  image de niveau l+1
 */
double** reduce(double** img, double** h, int size, int nl_old, int nc_old, int nl, int nc) {
    double** res = alloue_image_double(nl,nc);
    for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	    for (int ky=0; ky<size; ky++) {
    	        for (int kx=0; kx<size; kx++) {
    	            int y = 2*j+ky-(int)(size/2);
                    if (y<0) {
    	                y = nc_old+y;
    	            } else if (y>=nc_old) {
    	                y = y-nc_old;
                    }
                    int x = 2*i+kx-(int)(size/2);
                    if (x<0) {
    	                x = nl_old+x;
    	            } else if (x>=nl_old) {
    	                x = x-nl_old;
                    } 	      
    	            res[i][j] += img[x][y]*h[kx][ky];
                }
    	    }
    	}
    }
    return res;
}

/**
 * Expansion x4 d'une image de la pyramide de Gauss (par interpolation)
 */
double** expand(double** img, double** h, int size, int nl_old, int nc_old, int nl, int nc) {
    double** res = alloue_image_double(nl,nc);
    for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	    for (int n=0; n<size; n++) {
    	        for (int m=0; m<size; m++) {
    	            // Seulement si les indices sont des entiers
    	            if (((j-m+(int)(size/2)) % 2 == 0) && ((i-n+(int)(size/2)) % 2 == 0)) {
        	            int y = (j-m+(int)(size/2))/2;
                        if (y<0) {
        	                y = nc_old+y;
        	            } else if (y>=nc_old) {
        	                y = y-nc_old;
                        }
                        int x = (i-n+(int)(size/2))/2;
                        if (x<0) {
        	                x = nl_old+x;
        	            } else if (x>=nl_old) {
        	                x = x-nl_old;
                        } 	      
        	            res[i][j] += img[x][y]*h[n][m];
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
	int size = 9;
	double** h = getLowpassFilter();
	double** highpassFilter = getHighpassFilter();
    
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
    //double*** pyramidGaussExpand = = (double***) calloc(Nb,sizeof(double**));
    double*** pyramidLap = (double***) calloc(Nb,sizeof(double**));
    double*** reduceGauss = (double***) calloc(Nb+1,sizeof(double**));
    reduceGauss[0] = im1;
    
    double*** pyramidHighFilter = (double***) calloc(Nb,sizeof(double**));

    
    for (int k=0; k<Nb; k++) {
        reduceGauss[k+1] = reduce(reduceGauss[k], h, size, dimL[k], dimC[k], dimL[k+1], dimC[k+1]);
        pyramidGauss[k] = expand(reduceGauss[k+1], h, size, dimL[k+1], dimC[k+1], dimL[k], dimC[k]);
    }

    for (int k=0; k<Nb; k++) {
        pyramidLap[k] = alloue_image_double(dimL[k],dimC[k]);
        for (int i=0; i<dimL[k]; i++) {
    	    for (int j=0; j<dimC[k]; j++) {
    	        if (k==0) {
    	            pyramidLap[k][i][j] = im1[i][j]-pyramidGauss[0][i][j];
                } else {
                    pyramidLap[k][i][j] = reduceGauss[k][i][j]-pyramidGauss[k][i][j];
                }
            }
        }
    }


	for (int k=0; k<Nb; k++) {
    	if (k==0) {
    		pyramidHighFilter[0] = applyHighFilter(im1,highpassFilter,nl,nc);
    	} else {
    	    pyramidHighFilter[k] = applyHighFilter(reduceGauss[k],highpassFilter,dimL[k],dimC[k]);
    	}
    }

    
	// Ecriture
    for (int k=0; k<Nb; k++) {
        char* indice = (char*) calloc(1,sizeof(char));
        sprintf(indice, "%i", k+1);
        char filename1[20];
        strcpy(filename1, "Lowpass_level");
        strcat(filename1, indice);
        strcat(filename1, ".pgm");
        char filename2[20];
        strcpy(filename2, "Highpass_level");
        strcat(filename2, indice);
        strcat(filename2, ".pgm");
        char filename3[20];
        strcpy(filename3, "WithFilter_Highpass_level");
        strcat(filename3, indice);
        strcat(filename3, ".pgm");
        ecritureimagepgm(filename1, scale(pyramidGauss[k],dimL[k],dimC[k]), dimL[k], dimC[k]);
        ecritureimagepgm(filename2, scale(pyramidLap[k],dimL[k],dimC[k]), dimL[k], dimC[k]);
        ecritureimagepgm(filename3, scale(pyramidHighFilter[k],dimL[k],dimC[k]), dimL[k], dimC[k]);
    }
    
	for (int i=0; i<nl; i++) {
    	for (int j=0; j<nc; j++) {
    	printf("%f    %f\n", pyramidHighFilter[0][i][j], pyramidLap[0][i][j]);
    	
    }
    }
	
    return 0;
}  
