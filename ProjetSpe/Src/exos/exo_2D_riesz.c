// Exercice : retrouver Riesz en 2D en connaissant le résultat théorique

#include "pgm.h"
#include "fft.h"
#include <math.h>

int main (int ac, char **av) {

	double **image_input;
	double **image_input_imag;
	double **R1_theorique;
	double **R2_theorique;
	double **im2_reel;
    double **im2_imag;
    double **im3_reel;
    double **im3_imag;
    double **R1_reel;
    double **R1_imag;
    double **R2_reel;
    double **R2_imag;
    double **R1_reel_shift;
    double **R1_imag_shift;
    double **R2_reel_shift; 
    double **R2_imag_shift;
    double **R1_reel_spatial;
    double **R1_imag_spatial;
    double **R2_reel_spatial;
    double **R2_imag_spatial;
	
	// Parametres
	double beta = 1.5;
	double u0 = 10;
	
	// Nb de lignes et nb de colonnes (puissance de 2)
	int nl = 1024;
	int nc = 1024;
	
	double d = 5.0;
	double N = 1024.0;
	
	// Image input : f(x) = cos(u0(x1*cos(beta)+x2*sin(beta)))
	image_input = alloue_image_double(nl,nc);
    for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		image_input[u][v] = cos(u0*(d/N*(double)u*cos(beta)+d/N*(double)v*sin(beta)));
    	}              
    }
    unsigned char** affiche_input = scale(image_input,nl,nc);
    ecritureimagepgm(av[1], affiche_input, nl, nc);
    libere_image(affiche_input);
    
    // R1 et R2 théoriques
	R1_theorique = alloue_image_double(nl,nc);
	R2_theorique = alloue_image_double(nl,nc);
	for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		double sinus = sin(u0*(d/N*(double)u*cos(beta)+d/N*(double)v*sin(beta)));
    		R1_theorique[u][v] = cos(beta)*sinus;
    		R2_theorique[u][v] = sin(beta)*sinus;
    	}
    }
    unsigned char** affiche_R1_theorique = scale(R1_theorique,nl,nc);
    unsigned char** affiche_R2_theorique = scale(R2_theorique,nl,nc);
    ecritureimagepgm(av[2], affiche_R1_theorique, nl, nc);
    ecritureimagepgm(av[3], affiche_R2_theorique, nl, nc);
    libere_image(affiche_R1_theorique);
    libere_image(affiche_R2_theorique);
    
    
    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    
	// Experimentation
	
    image_input_imag = alloue_image_double(nl,nc);
    im2_reel = alloue_image_double(nl,nc);
    im2_imag = alloue_image_double(nl,nc);
    im3_reel = alloue_image_double(nl,nc);
    im3_imag = alloue_image_double(nl,nc);
    R1_reel = alloue_image_double(nl,nc);
    R1_imag = alloue_image_double(nl,nc);
    R2_reel = alloue_image_double(nl,nc);
    R2_imag = alloue_image_double(nl,nc);
    R1_reel_shift = alloue_image_double(nl,nc);
    R1_imag_shift = alloue_image_double(nl,nc);
    R2_reel_shift = alloue_image_double(nl,nc); 
    R2_imag_shift = alloue_image_double(nl,nc);
    R1_reel_spatial = alloue_image_double(nl,nc);
    R1_imag_spatial = alloue_image_double(nl,nc);
    R2_reel_spatial = alloue_image_double(nl,nc);
    R2_imag_spatial = alloue_image_double(nl,nc);
    
	
	// FFT de image_input -> im2
    fft(image_input, image_input_imag, im2_reel, im2_imag, nl, nc);
    
    // Translation : im2 -> im3
    fftshift(im2_reel, im2_imag, im3_reel, im3_imag, nl, nc);
    
    // Transformée de Riesz dans le domaine fréquentiel
    for (int u=0; u<nl; u++) {
    	double nu1 = ((double)u-N/2.0)/d;
    	for (int v=0; v<nc; v++) {
    		double nu2 = ((double)v-N/2.0)/d;
    		double norme = sqrt((double)(nu1*nu1+nu2*nu2));
    		if (norme != 0) {
    			R1_reel[u][v] = ((double) (-nu1*(double)im3_imag[u][v]) ) / norme;
    			R1_imag[u][v] = ((double) (+nu1*(double)im3_reel[u][v]) ) / norme;
    			R2_reel[u][v] = ((double) (-nu2*(double)im3_imag[u][v]) ) / norme;
    			R2_imag[u][v] = ((double) (+nu2*(double)im3_reel[u][v]) ) / norme;
    		}
    	}
    }
    
    // Translation inverse
    fftshift(R1_reel, R1_imag, R1_reel_shift, R1_imag_shift, nl, nc);
    fftshift(R2_reel, R2_imag, R2_reel_shift, R2_imag_shift, nl, nc);
    
    // IFFT
    ifft(R1_reel_shift, R1_imag_shift, R1_reel_spatial, R1_imag_spatial, nl, nc);
	ifft(R2_reel_shift, R2_imag_shift, R2_reel_spatial, R2_imag_spatial, nl, nc);
    
    ////////////////////////////////////////////////////////////////////////////


    unsigned char** affiche_R1_exp = scale(R1_reel_spatial,nl,nc);
    unsigned char** affiche_R2_exp = scale(R2_reel_spatial,nl,nc);
    ecritureimagepgm(av[4], affiche_R1_exp, nl, nc);
    ecritureimagepgm(av[5], affiche_R2_exp, nl, nc);
    libere_image(affiche_R1_exp);
    libere_image(affiche_R2_exp);
  
	// Libération des images
    libere_image_double(image_input);
	libere_image_double(image_input_imag);
	libere_image_double(im2_reel);
	libere_image_double(im2_imag);
	libere_image_double(im3_reel);
	libere_image_double(im3_imag);
	libere_image_double(R1_reel);
	libere_image_double(R1_imag);
	libere_image_double(R2_reel);
	libere_image_double(R2_imag);
	libere_image_double(R1_reel_shift);
	libere_image_double(R1_imag_shift);
	libere_image_double(R2_reel_shift);
	libere_image_double(R2_imag_shift);
	libere_image_double(R1_reel_spatial);
	libere_image_double(R1_imag_spatial);
	libere_image_double(R2_reel_spatial);
	libere_image_double(R2_imag_spatial);
	libere_image_double(R1_theorique);
	libere_image_double(R2_theorique);
    
    return 0;
}  
