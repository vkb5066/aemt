
#include <math.h>
#include <stdlib.h>
#include "FDecs.h"

//File for computing smearing parameters, interpolations

//Fills arr with a gaussian profile, updates its length
//Based on dx, the distance between two adjacent indices in ANY (not
//(just x) direction.  Cutoff = number of standard deviations before
//the profile is cut (3 to 4 is good)
void GaussianProfile(double **arr, int *arrLen, 
				     const double sigma, const int cutoff){
	int len = 2*(int)(ceil((double)cutoff*sigma) - 1.0);
	if(len%2 == 0){
		len++;
	}
	*arrLen = len;

	int mp = (len - 1)/2; ///the midpoint of the array
	double twoSigSqd = 2.0*sigma*sigma;
	*arr = malloc(len*sizeof(double));
	for(int i = 0; i < len; ++i){
		(*arr)[i] = 1.0/sqrt(PI*twoSigSqd) * 
			        exp(-(double)((i-mp)*(i-mp))/twoSigSqd);
	}
}

//Computes the cross-convolution between x and p(rofile), both of len
//len (should be an odd number), at x's midpoint
//Assumes that x and gProf are completly filled with appropriate 
//extrapolations, etc.  This is done in Grid.c
double CConvS(const double **x, const double *p, const int len){
	double ret = 0.0;
	for(unsigned int i = 0; i < len; ++i){ 
			ret += p[i]*(*x[i]);
	}

	return ret;
}


