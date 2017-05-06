#include "rieszPyramidFilters.h"

void affiche(double** h, int nl, int nc) {
	printf("Affichage d'un filtre\n");
	for (int i=0; i<nl; i++) {
		for (int j=0; j<nc; j++) {
			if (h[i][j] >= 0) {
				printf(" %f ", h[i][j]);
			} else {
				printf("%f ", h[i][j]);
			}
		}
		printf("\n");
	}
}


double** getLowpassFilter() {
	double** res = alloue_image_double(9,9);
	res[0][0] = -0.0001;
	res[1][0] = -0.0007;
	res[2][0] = -0.0023;
	res[3][0] = -0.0046;
	res[4][0] = -0.0057;
	res[1][1] = -0.0030;
	res[2][1] = -0.0047;
	res[3][1] = -0.0025;
	res[4][1] = -0.0003;
	res[2][2] = 0.0054;
	res[3][2] = 0.0272;
	res[4][2] = 0.0387;
	res[3][3] = 0.0706;
	res[4][3] = 0.0910;
	res[4][4] = 0.1138;
	for (int j=0; j<5; j++) {
		for (int i=4; i>j; i--) {
			res[j][i] = res[i][j];
		}
	}
	for (int j=5; j<9; j++) {
		for (int i=0; i<5; i++) {
			res[i][j] = res[i][8-j];
		}
	}
	for (int j=0; j<9; j++) {
		for (int i=5; i<9; i++) {
			res[i][j] = res[8-i][j];
		}
	}
	return res;
}


double** getHighpassFilter() {
	double** res = alloue_image_double(9,9);
	res[0][0] =0.0000;
	res[1][0] =0.0003;
	res[2][0] =0.0011;
	res[3][0] =0.0022;
	res[4][0] =0.0027;
	res[5][0] =0.0022;
	res[6][0] =0.0011;
	res[7][0] =0.0003;
	res[8][0] =0.0000;


	res[0][1] =0.0003 ;
	res[1][1] =0.0020 ;
	res[2][1] =0.0059 ;
	res[3][1] =0.0103  ;
	res[4][1] =0.0123 ;
	res[5][1] =0.0103 ;
	res[6][1] =0.0059 ;
	res[7][1] =0.0020 ;
	res[8][1] =0.0003 ;

	res[0][2] = 0.0011 ;
	res[1][2] =0.0059 ;
	res[2][2] = 0.0151 ;
	res[3][2] =0.0249 ;
	res[4][2] = 0.0292 ;
	res[5][2] = 0.0249 ;
	res[6][2] = 0.0151 ;
	res[7][2] = 0.0059 ;
	res[8][2] =0.0011 ;

	res[0][3] = 0.0022 ;
	res[1][3] =0.0103 ;
	res[2][3] = 0.0249 ;
	res[3][3] =0.0402 ;
	res[4][3] = 0.0469 ;
	res[5][3] = 0.0402 ;
	res[6][3] = 0.0249 ;
	res[7][3] = 0.0103 ;
	res[8][3] =0.0022 ;

	res[0][4] = 0.0027; 
	res[1][4] =0.0123; 
	res[2][4] = 0.0292; 
	res[3][4] =0.0469; 
	res[4][4] = -0.9455 ;
	res[5][4] = 0.0469; 
	res[6][4] = 0.0292 ;
	res[7][4] = 0.0123 ;
	res[8][4] =0.0027 ;

    for(int i=0; i<9; i++) {
   		for(int j =0;j<4;j++) {
	    	res[i][8-j] = res[i][j];
	    }
	}
	return res;
}
