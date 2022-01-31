//constants file

#include "FDecs.h"

//input-output variables
const extern int LINESIZE = 255;
const extern BASE = 10; //parameter to string-to-long int function

//Special values
const extern double E_LLIM = -1.0E50; ///lower limit of energy in eV
const extern double E_SPEC_DBL = -1.0E51; ///safe val for special 
										  ///values of double
const extern double EXP_MAX = +350.0; ///value used before 
									  ///exp(x) = inf
const extern double EXP_MIN = -350.0; ///value used before exp(x) = 0
const extern double ABS_MEFF_MIN = 1.0E-300;
const extern double ABS_MEFF_MAX = 1.0E300;

//Caluclation variables
const extern double PI = 3.141592654;
const extern double BOLTZMANN = 8.617333262145E-5; ///eV / K
const extern double MINV_REL_CONV = 0.131234202; ///angst^2/(eV s^2)
												 ///-> untiless
const extern double IDW_DIST_TOL = 1.0E-12; ///tol b4 gridpoint is 
						                    ///considered "on top of"
											///a k pt
const extern int STD_CUTOFF = 3;

