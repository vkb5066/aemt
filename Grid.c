

//Header file for the integration / derivative grid

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "FDecs.h"
#include "extern//kdtree-master//kdtree.h"

double Dist(const double ux, const double uy, const double uz, 
			const double vx, const double vy, const double vz){
	double diffX = vx-ux;
	double diffY = vy-uy;
	double diffZ = vz-uz;
	return sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
}

//Fills array arr with {dx, dy, dz}
void FillDifferentals(double **arr, const double *lims, const int *nGridPts){
	*arr = malloc(3*sizeof(double));
	for(unsigned int i = 0; i < 3; ++i){
		(*arr)[i] = (lims[2*i + 1] - lims[2*i])/((double)(nGridPts[i] - 1));
	}
}

//Gives the integration / derivative grid, as well as the connection
//tables / neighbor indices
struct gp** GiveInitGrid(const double *domain, const int *nGps,
						const void *rlKDNode, const double *kKDNode,
						const double kpRad, unsigned int *gridSize,
						int ***kCTable, double ***idTable, 
						int **kpdTableLens){
	///Find the index of the gamma point:
	///If any grid point's closest neighbor is gamma ({0., 0., 0.}),
	///then it is in the 1st brillouin zone
	int gIndex = GetGammaIndex(rlKDNode);

	///Setup tables
	unsigned int tableSize = nGps[0]*nGps[1]*nGps[2];

	double *allX = malloc(nGps[0]*sizeof(double));
	double *allY = malloc(nGps[1]*sizeof(double));
	double *allZ = malloc(nGps[2]*sizeof(double));
	Linspace(&allX, nGps[0], domain[0], domain[1]);
	Linspace(&allY, nGps[1], domain[2], domain[3]);
	Linspace(&allZ, nGps[2], domain[4], domain[5]);

	int **gCTable = malloc(tableSize*sizeof(int*));
	for(unsigned int i = 0; i < tableSize; i++){
		gCTable[i] = malloc(6*sizeof(int));
	}


	//Fill full rectangle around BZ - ignore the more intense calcs
	//until we have the reduced list (so just fill x, y, z, nbrs,
	//and determine if it is in bz)
	struct gp **fullList = malloc(tableSize*sizeof(struct gp*));
	int *indsToFree = calloc(tableSize, sizeof(int));

	///Internal variables (a waste to define these inside of a loop)
	///I'd be surprised if the compiler didnt take care of this, but
	///im doing it just in case
	int neighbors[6] = {-1, -1, -1, -1, -1, -1};
	int zOffsetN = -1; int zOffsetP = +1;
	int yOffsetN = -nGps[2]; int yOffsetP = +nGps[2];
	int xOffsetN = -nGps[2]*nGps[1]; int xOffsetP = +nGps[2]*nGps[1];
	double x = 0.; double y = 0.; double z = 0.;
	double dm[3] = {0., 0., 0.}; ///dummy variables so unix dosn't 
	void *lpRes = NULL;          ///freak out
	int lpIndex = 0;
	for(unsigned int i = 0, count = 0; i < nGps[0]; ++i){
	for(unsigned int j = 0; j < nGps[1]; ++j){
	for(unsigned int k = 0; k < nGps[2]; ++k, ++count){
		x = allX[i];
		y = allY[j];
		z = allZ[k];
				
		///grid point connection table
		gCTable[count][0] = count + xOffsetN;
		gCTable[count][1] = count + xOffsetP;
		gCTable[count][2] = count + yOffsetN;
		gCTable[count][3] = count + yOffsetP;
		gCTable[count][4] = count + zOffsetN;
		gCTable[count][5] = count + zOffsetP;
				

		///1st bz check
		lpRes = kd_nearest3(rlKDNode, x, y, z);
		lpIndex = (int)kd_res_item3(lpRes, &dm[0], &dm[1], &dm[2]);

		///Add this point to the list
		struct gp *pt = malloc(sizeof(struct gp));
		pt->crds[0] = x;
		pt->crds[1] = y;
		pt->crds[2] = z;
		pt->isInBz = (gIndex == lpIndex)? 't' : 'f'; 
		fullList[count] = pt;

		///Should this be in the loop?  Maybe declare voids 
		///outside of loop, update within, then free after
		kd_res_free(lpRes);
	}
	}
	}
	free(allX);
	free(allY);
	free(allZ);

	//Prep to complete connection tables
	*kCTable = calloc(tableSize, sizeof(int*));
	*idTable = calloc(tableSize, sizeof(double*));
	*kpdTableLens = calloc(tableSize, sizeof(int*));

	//(More) internal variables
	double x_ = 0.; double y_ = 0.; double z_ = 0.; ///k-point data
	char usePoint = '!';
	int testIndex = -1;

	//Extra bit to add to 1/(x) to avoid infinity
	const double extra = kpRad/(1.0E10);

	//Now reduce the list to points that are either in the BZ or 
	//close.  Do the complicated calculations for those in BZ
	struct gp **redArr = calloc(tableSize, sizeof(struct gp*));
	unsigned int count = 0;
	unsigned int freeCount = 0;
	for (unsigned int i = 0; i < tableSize; ++i){ 
		///If this point is in the BZ, keep it
		///If this point is not in the BZ, but has a neighbor that
		///is, keep it.  Otherwise, move on to the next point
		///Also, set up the neighbor pointers in the same loop
		struct gp *thisPt = fullList[i];
		usePoint = thisPt->isInBz;
		for(unsigned int j = 0; j < 6; ++j){
			testIndex = gCTable[i][j];

			if(testIndex < 0 || testIndex > tableSize - 1){
				thisPt->nbrs[j] = NULL;
				continue;
			}
			thisPt->nbrs[j] = fullList[testIndex];

			if (thisPt->nbrs[j]->isInBz == 't'){
				usePoint = 't';
			}
		}
		if(usePoint == 'f'){
			thisPt->isInBz = 'k'; ///k: that this will be freed later
			indsToFree[freeCount] = i; ///don't free now since this 
			freeCount++;               ///might be someone else's
			continue;                  ///neighbor
		}
		
		///If we get here, we're using this point
		x = fullList[i]->crds[0];
		y = fullList[i]->crds[1];
		z = fullList[i]->crds[2];


		///k-points: connection and distance tables
		struct kdres *kpRes = kd_nearest_range3(kKDNode, x, y, z, kpRad);
		int resSize = kd_res_size(kpRes);
		if (resSize == 0) {
			kpRes = kd_nearest3(kKDNode, x, y, z);
			resSize = 1;
		}

		(*kCTable)[count] = malloc(resSize*sizeof(int));
		(*idTable)[count] = malloc(resSize*sizeof(double));
		(*kpdTableLens)[count] = resSize;

		for (unsigned int j = 0; j < resSize; ++j) {
			////index + update coords of this k point
			int thisKIndex = (int*)kd_res_item3(kpRes, &x_, &y_, &z_);

			(*kCTable)[count][j] = thisKIndex;
			(*idTable)[count][j] = 1.0/(Dist(x, y, z, x_, y_, z_) + 
										extra);

			kd_res_next(kpRes);
		}

		//Finish initializing, clean up
		redArr[count] = thisPt;
		kd_res_free(kpRes);
		count++;
	}

	//Clean up
	for(unsigned int i = 0; i < tableSize; ++i){
		if(i > count){
			free(redArr[i]);
			free((*kCTable)[i]);
			free((*idTable)[i]);
			free((*kpdTableLens)[i]);
		}
		free(gCTable[i]);
	}
	free(gCTable);
	for(unsigned int i = 0; i < count; ++i){
		for(unsigned int j = 0; j < 6; ++j){
			if(redArr[i]->nbrs[j]->isInBz == 'k'){
				redArr[i]->nbrs[j] = NULL;
			}
		}
	}
	for(unsigned int i = 0; i < freeCount; ++i){
		free(fullList[indsToFree[i]]);
	}
	free(indsToFree);
	free(fullList);

	*gridSize = count;
	return redArr;
}

//Set energies using inverse distance weighting
//Meant to be multithreaded - beg and end are the beginning and end
//indices of the grid for this thread to go over.  There are no 
//mutexs here, so make sure begs, ends never overlap!!!
void SetEnergies(struct gp ***grid, const int beg, const int end,
				 const int idwPow, const double rad, 
				 const int bandIndex, 
				 const int *tableLens, const double **iDists, 
				 const double **eigs, const int **conTable){

	double num = 0.0; double den = 0.0;
	double wght = 0.0;
	const double iRad = 1.0/rad;
	for(unsigned int i = beg; i < end; ++i){
		///Inverse distance weighting: doi:10.1145/800186.810616
		///Update 2.7.2022: Improved for small choice of idwRad
		num = 0.0; den = 0.0;
		for(unsigned int j = 0; j < tableLens[i]; ++j){
			if(iDists[i][j] > INV_IDW_DIST_TOL){
				num = eigs[conTable[i][j]][bandIndex];
				den = 1.0;
				break;
			}

			wght = IntPow_d(iDists[i][j] - iRad, idwPow);
			num += wght*eigs[conTable[i][j]][bandIndex];
			den += wght;
		}

		(*grid)[i]->nrg = num/den;
	}
}

//Sets the first derivative for all grid[beg] to grid[end - 1]
//dir is the direction to take: 0, 1, 2 for x, y, z
void SetSmoothEnergies(struct gp ***grid, const int beg, 
					   const int end, const double *sigmas, 
					   const double *diffs){
	for(int dir = 0; dir < 3; ++dir){
		///Sigma small: no gaussian smoothing case
		if(sigmas[dir] < 1E-5){
			for(unsigned int i = beg; i < end; ++i){
				(*grid)[i]->nrgS[dir] = (*grid)[i]->nrg;
			}
			continue;
		}

		int lo = 2*dir; int hi = 2*dir + 1;
	
		///Gaussian profile, cross-conv array setup
		double *gProf; int gProfLen;
		GaussianProfile(&gProf, &gProfLen, sigmas[dir], STD_CUTOFF);
		int midpt = (gProfLen - 1)/2;

		double **x = malloc(gProfLen*sizeof(double*));
		double **xOrig = malloc(gProfLen*sizeof(double*));
		double **toFree = malloc(2*gProfLen*sizeof(double*));
		for(unsigned int i = 0; i < gProfLen; ++i){
			x[i] = malloc(sizeof(double));
			xOrig[i] = malloc(sizeof(double));
			toFree[2*i] = x[i];
			toFree[2*i + 1] = xOrig[i];
		}
	
		for(unsigned int i = beg; i < end; ++i){
			///fill x with k point indices from low to high, and
			///grid[i] right in the middle to cross-convolute
			memcpy(x, xOrig, gProfLen*sizeof(double*));

			*x[midpt] = (*grid)[i]->nrg; ////midpoint
			struct gp* curr = (*grid)[i];
			////low
			for(int j = midpt - 1, ext = 1; j > -1; --j){
				if(curr->nbrs[lo] != NULL){
					curr = curr->nbrs[lo];
					*x[j] = curr->nrg;
					continue;
				}

				x[j] = x[j + ext]; //mirror about BZ edge
				ext += 2;
			}
			curr = (*grid)[i];
			////high
			for(int j = midpt + 1, ext = 1; j < gProfLen; ++j){
				if(curr->nbrs[hi] != NULL){
					curr = curr->nbrs[hi];
					*x[j] = curr->nrg;
					continue;
				}
				
				x[j] = x[j - ext]; //mirror about BZ edge
				ext += 2;
			}
		
			///Set smoothed energy
			(*grid)[i]->nrgS[dir] = CConvS(x, gProf, gProfLen);
		}

		free(gProf);
		for(unsigned int i = 0; i < 2*gProfLen; ++i){
			free(toFree[i]);
		}
		free(x);
		free(xOrig);
		free(toFree);
	}

	//Set the absolute energies to average of the smoothed energies
	//since they'll be used for fermi distrubutions.  This introduces
	//a small error in the E - mu term (the band extrema shift lower/
	//higher for valance / conduction bands.
	//
	// 	   TODO:
	// 
	//Attempt to fix this by keeping track of the original extrema 
	//value then adding a constant correction term to each band's
	//smoothed energy.  This is easy if [beg, end) contains the whole
	//grid (as it does for no multithreading), but if multithreading
	//is implemented, [beg, end) may not contain the band extrema.  
	for(unsigned int i = beg; i < end; ++i){
		(*grid)[i]->nrg = ((*grid)[i]->nrgS[0] + 
						   (*grid)[i]->nrgS[1] + 
						   (*grid)[i]->nrgS[2])/3.0;
	}
}

//Sets the first derivatives in the direction indicated by dir by a 
//second order polynomial fit (f/b diffs are unstable by comparison)
//0, 1, 2 = x, y, z
//only sets derivatives for points in the BZ
void SetDerivs(struct gp*** grid, const int beg, const int end){
	for(int dir = 0; dir < 3; ++dir){
		int lo = 2*dir; int hi = 2*dir + 1;

		double a = 0.; double b = 0.;
		double xLo = 0.; double xMi = 0.; double xHi = 0.;
		double yLo = 0.; double yMi = 0.; double yHi = 0.;
		for(unsigned int i = beg; i < end; ++i){
			if((*grid)[i]->isInBz != 't'){
				continue;
			}

			///If a point is in the BZ, its "guarenteed" that it has
			///lo, hi neighbors - "no need to check for NULL"
			///**come back here when this proves to be wrong**
			xLo = (*grid)[i]->nbrs[lo]->crds[dir];
			yLo = (*grid)[i]->nbrs[lo]->nrgS[dir];
			xMi = (*grid)[i]->crds[dir];
			yMi = (*grid)[i]->nrgS[dir];
			xHi = (*grid)[i]->nbrs[hi]->crds[dir];
			yHi = (*grid)[i]->nbrs[hi]->nrgS[dir];
			PolyAB(&a, &b, xLo, yLo, xMi, yMi, xHi, yHi);

			(*grid)[i]->nrgDeriv[dir] = 2.0*a*xMi + b;
		}
	}
}

//Integrates over grid's BZ and fills results as:
//results = {xx, yy, zz, yz, zx, xy, fermi}
//if any = 0, that means there may have been numerical problems
//(i.e. exp^x was too big, etc...)
//This is for the conduction band, i.e. electrons
void IntegrateC(const struct gp **grid, const int beg, 
				const int end, double **results, const double t, 
				const double mu, const double *diffs){
	*results = calloc(7, sizeof(double));

	double dV = diffs[0]*diffs[1]*diffs[2];
	double kt = BOLTZMANN*t;
	for(unsigned int i = beg; i < end; ++i){
		if(grid[i]->isInBz != 't'){
			continue;
		}
		double dFermiDeriv = FDDeriv(grid[i]->nrg, mu, kt)*dV;

		///numerators
		(*results)[0] -= grid[i]->nrgDeriv[0]*grid[i]->nrgDeriv[0]*
					     dFermiDeriv;
		(*results)[1] -= grid[i]->nrgDeriv[1]*grid[i]->nrgDeriv[1]*
						 dFermiDeriv;
		(*results)[2] -= grid[i]->nrgDeriv[2]*grid[i]->nrgDeriv[2]*
						 dFermiDeriv;
		(*results)[3] -= grid[i]->nrgDeriv[1]*grid[i]->nrgDeriv[2]*
						 dFermiDeriv;
		(*results)[4] -= grid[i]->nrgDeriv[2]*grid[i]->nrgDeriv[0]*
						 dFermiDeriv;
		(*results)[5] -= grid[i]->nrgDeriv[0]*grid[i]->nrgDeriv[1]*
						 dFermiDeriv;

		///denominator
		(*results)[6] += FD(grid[i]->nrg, mu, kt)*dV;
	}
}

//Integrates over grid's BZ and fills results as:
//results = {xx, yy, zz, yz, zx, xy, fermi}
//if any = 0, that means there may have been numerical problems
//(i.e. exp^x was too big, etc...)
//This is for the valance band, i.e. holes
void IntegrateV(const struct gp** grid, const int beg, 
				const int end, double** results, const double t, 
				const double mu, const double* diffs){
	*results = calloc(7, sizeof(double));

	double dV = diffs[0]*diffs[1]*diffs[2];
	double kt = BOLTZMANN*t;
	for(unsigned int i = beg; i < end; ++i){
		if(grid[i]->isInBz != 't'){
			continue;
		}
		double dFermiDeriv = -FDDeriv(grid[i]->nrg, mu, kt)*dV;

		///numerators
		(*results)[0] += grid[i]->nrgDeriv[0]*grid[i]->nrgDeriv[0]*
						 dFermiDeriv;
		(*results)[1] += grid[i]->nrgDeriv[1]*grid[i]->nrgDeriv[1]*
						 dFermiDeriv;
		(*results)[2] += grid[i]->nrgDeriv[2]*grid[i]->nrgDeriv[2]*
						 dFermiDeriv;
		(*results)[3] += grid[i]->nrgDeriv[1]*grid[i]->nrgDeriv[2]*
						 dFermiDeriv;
		(*results)[4] += grid[i]->nrgDeriv[2]*grid[i]->nrgDeriv[0]*
						 dFermiDeriv;
		(*results)[5] += grid[i]->nrgDeriv[0]*grid[i]->nrgDeriv[1]*
						 dFermiDeriv;

		///denominator
		(*results)[6] += (1.0 - FD(grid[i]->nrg, mu, kt))*dV;
	}
}