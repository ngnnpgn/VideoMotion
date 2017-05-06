// Exercice : gradient en 2D

#include "pgm.h"
#include "fft.h"
#include <math.h>

#define PI 3.14159265358979323846

int main (int ac, char **av) {

	double **image_input;
	double **image_input_imag;
	double **Gx_theorique;
	double **Gy_theorique;
	double **im2_reel;
    double **im2_imag;
    double **im3_reel;
    double **im3_imag;
    double **Gx_reel;
    double **Gx_imag;
    double **Gy_reel;
    double **Gy_imag;
    double **Gx_reel_shift;
    double **Gx_imag_shift;
    double **Gy_reel_shift; 
    double **Gy_imag_shift;
    double **Gx_reel_spatial;
    double **Gx_imag_spatial;
    double **Gy_reel_spatial;
    double **Gy_imag_spatial;

	
	// Nb de lignes et nb de colonnes (puissance de 2)
	int nl = 512;
	int nc = 512;
	
	double d = 5.0;
	double N = 512.0;

	// Image input : f(x) = cos(2*pi*x+2*pi*y)
	image_input = alloue_image_double(nl,nc);
    for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		image_input[u][v] = cos(2*PI*(d/N*(double)u)+2*PI*(d/N*(double)v));
    	}
    }

    // Gx et Gy théoriques
	Gx_theorique = alloue_image_double(nl,nc);
	Gy_theorique = alloue_image_double(nl,nc);
	for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		Gx_theorique[u][v] = -2.0*PI*sin(2*PI*(d/N*(double)u)+2*PI*(d/N*(double)v));
    		Gy_theorique[u][v] = -2.0*PI*sin(2*PI*(d/N*(double)u)+2*PI*(d/N*(double)v));
    	}
    }
    
    /**************************************************************************/
    
	// Experimentation
	
    image_input_imag = alloue_image_double(nl,nc);
    im2_reel = alloue_image_double(nl,nc);
    im2_imag = alloue_image_double(nl,nc);
    im3_reel = alloue_image_double(nl,nc);
    im3_imag = alloue_image_double(nl,nc);
    Gx_reel = alloue_image_double(nl,nc);
    Gx_imag = alloue_image_double(nl,nc);
    Gy_reel = alloue_image_double(nl,nc);
    Gy_imag = alloue_image_double(nl,nc);
    Gx_reel_shift = alloue_image_double(nl,nc);
    Gx_imag_shift = alloue_image_double(nl,nc);
    Gy_reel_shift = alloue_image_double(nl,nc); 
    Gy_imag_shift = alloue_image_double(nl,nc);
    Gx_reel_spatial = alloue_image_double(nl,nc);
    Gx_imag_spatial = alloue_image_double(nl,nc);
    Gy_reel_spatial = alloue_image_double(nl,nc);
    Gy_imag_spatial = alloue_image_double(nl,nc);
    
	// FFT de image_input -> im2
    fft(image_input, image_input_imag, im2_reel, im2_imag, nl, nc);
    
    // Translation : im2 -> im3
    fftshift(im2_reel, im2_imag, im3_reel, im3_imag, nl, nc);
    
    // Transformée de Riesz dans le domaine fréquentiel
    for (int u=0; u<nl; u++) {
    	double nu1 = ((double)u-N/2.0)/d;
    	for (int v=0; v<nc; v++) {
    		double nu2 = ((double)v-N/2.0)/d;
    		Gx_reel[u][v] = (double) ( 2*PI*nu1*(double)im3_imag[u][v]);
    		Gx_imag[u][v] = (double) (-2*PI*nu1*(double)im3_reel[u][v]);
    		Gy_reel[u][v] = (double) ( 2*PI*nu2*(double)im3_imag[u][v]);
    		Gy_imag[u][v] = (double) (-2*PI*nu2*(double)im3_reel[u][v]);
    	}
    }
    
    // Translation inverse
    fftshift(Gx_reel, Gx_imag, Gx_reel_shift, Gx_imag_shift, nl, nc);
    fftshift(Gy_reel, Gy_imag, Gy_reel_shift, Gy_imag_shift, nl, nc);
    
    // IFFT
    ifft(Gx_reel_shift, Gx_imag_shift, Gx_reel_spatial, Gy_imag_spatial, nl, nc);
	ifft(Gy_reel_shift, Gy_imag_shift, Gy_reel_spatial, Gy_imag_spatial, nl, nc);
    
    /**************************************************************************/

    // Comparaison
    
    for (int u=0; u<nl; u++) {
    	for (int v=0; v<nc; v++) {
    		if (fabs(Gx_reel_spatial[u][v]-Gx_theorique[u][v])>0.03) {
    		    printf("ERREUR\n");
    		}
    		if (fabs(Gy_reel_spatial[u][v]-Gy_theorique[u][v])>0.03) {
    		    printf("ERREUR\n");
    		}
    	}
    } 
  
    printf("FIN\n");
    
	// Libération des images
	libere_image_double(image_input);
	libere_image_double(image_input_imag);
	libere_image_double(im2_reel);
	libere_image_double(im2_imag);
	libere_image_double(im3_reel);
	libere_image_double(im3_imag);
	libere_image_double(Gx_reel);
	libere_image_double(Gx_imag);
	libere_image_double(Gy_reel);
	libere_image_double(Gy_imag);
	libere_image_double(Gx_reel_shift);
	libere_image_double(Gx_imag_shift);
	libere_image_double(Gy_reel_shift);
	libere_image_double(Gy_imag_shift);
	libere_image_double(Gx_reel_spatial);
	libere_image_double(Gx_imag_spatial);
	libere_image_double(Gy_reel_spatial);
	libere_image_double(Gy_imag_spatial);
	libere_image_double(Gx_theorique);
	libere_image_double(Gy_theorique);
    
    return 0;
}  
