

// f(x) = cos ( x1 + x2 )

#include "pgm.h"
#include <math.h>

int main (int ac, char **av) {
	
	double beta = 0.6;
	double u0 = 0.07;
	
    // Variables images
	unsigned char **res=NULL;
	double **input=NULL;
	double **R1=NULL;
	double **R2=NULL;
	double **A=NULL;

    // Nb de lignes et nb de colonnes
	int nl = 512;
	int nc = 512;

	// Vérification du nombre de paramètres
    if (ac < 4) {
        printf("Usage : %s input R1 R2\n", av[0]); exit(1);
    }
    
	res = alloue_image(nl,nc);
	input = alloue_image_double(nl,nc);
	R1 = alloue_image_double(nl,nc);
	R2 = alloue_image_double(nl,nc);
	A = alloue_image_double(nl,nc);
	
	/* Input */
    for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		input[u][v] = (double) cos( u0*((double)u*cos(beta) + (double)v*sin(beta)) );
    	//	res[u][v] = (input[u][v]+1.0)/2.0*255.0;
    	}
    }
    double **inter = scale(input,nl,nc);
    res = imdouble2uchar(inter,nl,nc);
    ecritureimagepgm(av[1], res, nl, nc);

	/* R1 */
	for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		R1[u][v] = (double) cos(beta)*sin( u0*((double)u*cos(beta) + (double)v*sin(beta)) );
    		res[u][v] = (R1[u][v]+1.0)/2.0*255.0;
    	}
    }
    ecritureimagepgm(av[2], res, nl, nc);

	/* R2 */
	for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		R2[u][v] = (double) sin(beta)*sin( u0*((double)u*cos(beta) + (double)v*sin(beta)) );
    		res[u][v] = (R2[u][v]+1.0)/2.0*255.0;
    	}
    }
    ecritureimagepgm(av[3], res, nl, nc);
    
    /* Orientation */
    for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		double arctan = atan2(sqrt(R2[u][v]*R2[u][v]+R1[u][v]*R1[u][v]),input[u][v]);
    		res[u][v] = (arctan+3.14)/(2.0*3.14) * 255.0;
    	}
    }
    ecritureimagepgm(av[4], res, nl, nc);
    
    
    /* Amplitude */
    for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		A[u][v] = sqrt(R2[u][v]*R2[u][v]+R1[u][v]*R1[u][v]+input[u][v]*input[u][v]);
    	}
    }
    //res = imdouble2uchar(scale(A,nl,nc),nl,nc);
    ecritureimagepgm(av[5], imdouble2uchar(A,nl,nc), nl, nc);
   

	// Libération mémoire
    libere_image(res); 
    
    return 0;
}  
