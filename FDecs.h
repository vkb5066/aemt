//Function declarations for AEM code

//External variables
const extern int LINESIZE;
const extern int BASE;

const extern double E_LLIM;
const extern double E_SPEC_DBL;
const extern double EXP_MAX;
const extern double EXP_MIN;
const extern double ABS_MEFF_MIN;
const extern double ABS_MEFF_MAX;

const extern double PI;
const extern double BOLTZMANN;
const extern double MINV_REL_CONV;
const extern double INV_IDW_DIST_TOL; 
const extern int STD_CUTOFF;

//Input / output
struct minv{
	double numer[6]; //sum(integral(vi*vj*dFd/dE*dV))
	double denom; //sum(integral(Fd*dV))
	double temp;
	double mu;
	char type;
};
void FillArray_d(double**, const int, const double, const double);
void FillArray_i(int**, const int, const int, const int);

void ReadInputFile(const char*, double**, int**, double***, int*, 
				   double*, double**, unsigned int*, double**, 
				   unsigned int*, double**, unsigned int*, double**, 
				   unsigned int*, double**, unsigned int*, int**, 
				   unsigned int*, int**);
void ReadEigennvalueFile(const char*, int*, double***, int*, 
						 double***);
void WriteResults(const char*, const struct minv***, const int, 
				  const int);

//Grid
struct gp{
	double crds[3]; ///{x, y, z}
	double nrg; ///raw energy (for fermi distrib)
	double nrgS[3]; ///smoothed energy (for derivs): {x, y, z}
	double nrgDeriv[3]; ///{x, y, z}
	struct gp *nbrs[6]; ///{-x, +x, -y, +y, -z, +z}
	char isInBz;
};
double Dist(const double, const double, const double, const double, 
			const double, const double);
void FillDifferentals(double**, const double*, const int*);
struct gp** GiveInitGrid(const double*, const int*, const void*, 
						 const double*, const double, unsigned int*, 
						 int***, double***, int**);
void SetEnergies(struct gp***, const int, const int, const int, 
				 const double, const int, const int*, const double**, 
				 const double**, const int**);
void SetSmoothEnergies(struct gp***, const int, const int, 
					   const double*, const double*);
void SetDerivs(struct gp***, const int, const int);
void IntegrateC(const struct gp**, const int, const int, double**, 
			    const double, const double, const double*);
void IntegrateV(const struct gp**, const int, const int, double**,
				const double, const double, const double*);


//Gaussian Smearing
void GaussianProfile(double**, int*, const double, 
					 const int);
double CConvS(const double**, const double*, const int);


//Helper functions
void SetKdTree(void**, const int, const double**);
void FreeArr_nxm_i(int***, const int);
void FreeArr_nxm_d(double***, const int);
void Linspace(double**, const int, const double, const double);
void SplitArr(const struct gp**, const int, const int, 
			  unsigned int***);
void InitMinv(struct minv****, const int, const double*, const int, 
			  const double*);
void FillMinv(struct minv****, const double*, const char, const int, 
			  const int);
void FreeMinv(struct minv****, const int, const int);
void ClearStdout();

//Math
int IntPow_i(const int, const int);
double IntPow_d(const double, const int);
void PolyAB(double*, double*, const double, const double, 
		    const double, const double, const double, 
			const double);
double FD(const double, const double, const double);
double FDDeriv(const double, const double, const double);
double MInvConv(const double, const double);

//Recip lattice
double** GetRecipLatPts(const int, const double**);
int GetGammaIndex(const void*);

//Extra tools
double DistSqd(const double, const double, const double,
			   const struct gp*);
double DistArr(const double*, const double*);
int FindIndex(const struct gp**, const int, const double, 
			  const double, const double);
int FindMaxima(const struct gp**, const int, const int);
int FindMinima(const struct gp**, const int, const int);
void WriteBS(const char*, const struct gp**, const int, const int);
void WriteDiagBS(const char*, const struct gp**, const int, 
				 const int, const int);
void WriteFullBs(const char*, const struct gp**, const unsigned int);