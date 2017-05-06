/* Premier essai : transformée de Riesz par FFT */

#include "pgm.h"
#include "fft.h"
#include <math.h>

#define PI 3.14159265358979323846

int main (int ac, char **av) {

    // Variables images
	unsigned char **im1=NULL, **res1=NULL, **res2=NULL, **res1_bis=NULL, **res2_bis=NULL, **result;
	double **im1_double;
	double **im1_reel, **im1_imag;
	double **im2_reel, **im2_imag;
	double **im3_reel, **im3_imag;
	double **R1_reel, **R1_imag;
	double **R2_reel, **R2_imag;
	double **R1_reel_shift, **R1_imag_shift;
	double **R2_reel_shift, **R2_imag_shift;
	double **R1_reel_spatial, **R1_imag_spatial;
	double **R2_reel_spatial, **R2_imag_spatial;
	// En coordonnées sphériques
	double **orientation;
	double **amplitude;
	double **phase;
	unsigned char **ori, **ori_bis;
	unsigned char **ampli, **ampli_bis;
	unsigned char **phi, **phi_bis;

    // Nb de lignes et nb de colonnes
	int nl,nc;
	// Nb de lignes et nb de colonnes en puissance de 2
	int nl2, nc2;

	// Vérification du nombre de paramètres
    if (ac < 4) {
        printf("Usage : %s entree R1 R2 Amplitude Orientation Phase\n", av[0]); exit(1);
    }
    
    // Lecture d'une image pgm dont le nom est passé sur la ligne de commande
    im1 = lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL) {
        puts("Lecture image impossible");
        exit(1);
    }
    nl2 = nl;
    nc2 = nc;
    
	// Conversion de l'image de unsigned char** en double**
    im1_double = imuchar2double(im1,nl,nc);

    // Pour que nl et nc soient des puissances de 2
    im1_reel = padimdforfft(im1_double,&nl2,&nc2);

    // Allocation des images (nulles)
    im1_imag = alloue_image_double(nl2,nc2); 
    im2_reel = alloue_image_double(nl2,nc2); im2_imag = alloue_image_double(nl2,nc2);
    im3_reel = alloue_image_double(nl2,nc2); im3_imag = alloue_image_double(nl2,nc2);
    
    R1_reel = alloue_image_double(nl2,nc2); R1_imag = alloue_image_double(nl2,nc2);
    R2_reel = alloue_image_double(nl2,nc2); R2_imag = alloue_image_double(nl2,nc2);
    
    R1_reel_shift = alloue_image_double(nl2,nc2); R1_imag_shift = alloue_image_double(nl2,nc2);
    R2_reel_shift = alloue_image_double(nl2,nc2); R2_imag_shift = alloue_image_double(nl2,nc2);
    
    R1_reel_spatial = alloue_image_double(nl2,nc2);
    R1_imag_spatial = alloue_image_double(nl2,nc2);
    
    R2_reel_spatial = alloue_image_double(nl2,nc2);
    R2_imag_spatial = alloue_image_double(nl2,nc2);
    
    orientation = alloue_image_double(nl2,nc2);
	ori = alloue_image(nl2,nc2);
	ori_bis = alloue_image(nl,nc);
	
	amplitude = alloue_image_double(nl2,nc2);
	ampli = alloue_image(nl2,nc2);
	ampli_bis = alloue_image(nl,nc);
	
	phase = alloue_image_double(nl2,nc2);
	phi = alloue_image(nl2,nc2);
	phi_bis = alloue_image(nl,nc);
	
	/**************************************************************************/
	/**************************************************************************/
	/**************************************************************************/
	
	// FFT de im1 -> im2
    fft(im1_reel, im1_imag, im2_reel, im2_imag, nl2, nc2);
    
    // Translation : im2 -> im3
    fftshift(im2_reel, im2_imag, im3_reel, im3_imag, nl2, nc2);
    
    double nld = (double) nl2;
    double ncd = (double) nc2;
    // Transformée de Riesz dans le domaine fréquentiel
    for (int u=0; u<nl2; u++) {
    	double nu1 = ((double)u-nld/2.0)/nld;
    	for (int v=0; v<nc2; v++) {
    		double nu2 = ((double)v-ncd/2.0)/ncd;
    		double norme = sqrt( (double) (nu1*nu1+nu2*nu2) );
    		if ( norme == 0) {
    			R1_reel[u][v] = -im3_imag[u][v];
    			R1_imag[u][v] = im3_reel[u][v];
    			R2_reel[u][v] = -im3_imag[u][v];
    			R2_imag[u][v] = im3_reel[u][v];
    		} else {

    			R1_reel[u][v] = ((double) (-nu1*(double)im3_imag[u][v]) ) / norme;
    			R1_imag[u][v] = ((double) (+nu1*(double)im3_reel[u][v]) ) / norme;
    			R2_reel[u][v] = ((double) (-nu2*(double)im3_imag[u][v]) ) / norme;
    			R2_imag[u][v] = ((double) (+nu2*(double)im3_reel[u][v]) ) / norme;

    		}
    	}
    }
    
    // Translation
    fftshift(R1_reel, R1_imag, R1_reel_shift, R1_imag_shift, nl2, nc2);
    fftshift(R2_reel, R2_imag, R2_reel_shift, R2_imag_shift, nl2, nc2);
    
    // IFFT
    ifft(R1_reel_shift, R1_imag_shift, R1_reel_spatial, R1_imag_spatial, nl2, nc2);
    ifft(R2_reel_shift, R2_imag_shift, R2_reel_spatial, R2_imag_spatial, nl2, nc2);
    
    ////////////////////////////////////////////////////////////////////////////

    /* Orientation */
    for (int u=0; u<nl2; u++) {
    	for (int v=0; v<nc2; v++) {
    		orientation[u][v] = atan2(R2_reel_spatial[u][v],R1_reel_spatial[u][v]);
    	}
    }

    res1 = scale(R1_reel_spatial, nl2, nc2);
    res2 = scale(R2_reel_spatial, nl2, nc2);
    ori = scale(orientation,nl2,nc2);

    /* Amplitude */
    for (int u=0; u<nl2; u++) {
    	for (int v=0; v<nc2; v++) {
    	    double r1 = (double) res1[u][v];
    	    double r2 = (double) res2[u][v];
    	    double in = (double) im1_reel[u][v];
    	    double A = sqrt(r1*r1+r2*r2+in*in);
    		amplitude[u][v] = A;
    		phase[u][v] = atan2(A,in);
    	}
    }
    ampli = scale(amplitude,nl2,nc2);
    phi = scale(phase,nl2,nc2);
        
	// Retour à nl et nc
    res1_bis = crop(res1,0,0,nl,nc);
    res2_bis = crop(res2,0,0,nl,nc);
    ori_bis =  crop(ori,0,0,nl,nc);
    ampli_bis = crop(ampli,0,0,nl,nc);
    phi_bis = crop(phi,0,0,nl,nc);

    // Ecritures des images dans les fichiers de sortie
        ecritureimagepgm(av[2], res1_bis, nl, nc);
	ecritureimagepgm(av[3], res2_bis, nl, nc);
	ecritureimagepgm(av[4], ori_bis, nl, nc);
	ecritureimagepgm(av[5], ampli_bis, nl, nc);
	ecritureimagepgm(av[6], phi_bis, nl, nc); 
	
	// Libération des images
	libere_image_double(im1_double);
	libere_image_double(im1_imag);
	libere_image_double(im2_reel);
	libere_image_double(im2_imag);
	libere_image_double(R1_reel);
	libere_image_double(R1_imag);
	libere_image_double(R1_reel_spatial);
	libere_image_double(R1_imag_spatial);
	libere_image_double(R2_reel);
	libere_image_double(R2_imag);
	libere_image_double(R2_reel_spatial);
	libere_image_double(R2_imag_spatial);
    libere_image(im1);
    libere_image(res1);
    libere_image(res2);
    libere_image(res1_bis);
    libere_image(res2_bis);
    
    return 0;
}  
