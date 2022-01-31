//Header for general math
//i.e. less safe but faster implementations of some common functions

#include <math.h>
#include "FDecs.h"

//Raise an int to an integer power
int IntPow_i(const int val, const int pow){
	int sum = val;
	for(unsigned int i = 1; i < pow; ++i){
		sum *= val;
	}
	return sum;
}
//Raise a double to an integer power
double IntPow_d(const double val, const int pow) {
	double sum = val;
	for (unsigned int i = 1; i < pow; ++i) {
		sum *= val;
	}
	return sum;
}

//Fits a second order polynomial to three points and returns the
//a and b values (c is unnecessary for derivs, and isn't super 
//accurate since the y's will be gaussian-smoothed)
//*1 is low, *2 is mid, and *3 is high
void PolyAB(double *a, double *b,
			const double x1, const double y1, const double x2, 
			const double y2, const double x3, const double y3){
	*a=(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))/((x1-x2)*(x1-x3)*(x2-x3));
	*b=(y2-y1)/(x2-x1)-(*a)*(x1+x2);
}

//Fermi dirac distribution for an energy, chemical potential, and a
//thermal energy
double FD(const double nrg, const double mu, const double kt){
	double x = (nrg - mu)/kt;
	if(x < EXP_MIN){
		return 1.0;
	}
	if(x > EXP_MAX){
		return 0.0;
	}

	return 1.0/(exp(x) + 1.0);
}

//Fermi dirac distribution deriv w.r.t. energy
double FDDeriv(const double nrg, const double mu, const double kt){
	double x = (nrg - mu)/kt;
	if (x < EXP_MIN) {
		return 0.0;
	}
	if (x > EXP_MAX) {
		return 0.0;
	}

	double expXp1 = exp(x) + 1.0;
	return -1.0/kt * (expXp1 - 1.0)/(expXp1*expXp1);
}

//Converts inverse effective mass from angst^2/(eV s^2) -> unitless
//(i.e. convert to 1/kg then multiply by m_0)
double MInvConv(const double num, const double den){
	double aNum = fabs(num); double aDen = fabs(den);
	if(aNum < ABS_MEFF_MIN || aDen < ABS_MEFF_MIN){
		return -0.00;
	}
	if (aNum > ABS_MEFF_MAX || aDen > ABS_MEFF_MAX){
		return -0.00;
	}

	return aNum/aDen*MINV_REL_CONV;
}