#pragma warning(disable:4996) //fuck you - why would throw an ERROR 
                              //for an unsafe standard io function?

//Debugging tools, etc.  Not to be used in the (final)
//main function
#include <stdio.h>
#include <stdlib.h>
#include "FDecs.h"

const int MAX_SIZE = 1024;

//Helper function to FindIndex(): distance between point at x, y, z,
//and a grid point
double DistSqd(const double x, const double y, const double z,
			   const struct gp *gridpt){
	double dx = x - gridpt->crds[0];
	double dy = y - gridpt->crds[1];
	double dz = z - gridpt->crds[2];
	return dx*dx + dy*dy + dz*dz;
}

//Finds the closest grid point to x, y, z and returns its index
int FindIndex(const struct gp **grid, const int gridSize, 
			  const double x, const double y, const double z){
	int bestInd = 0;
	double bestDist = DistSqd(x, y, z, grid[0]);
	for(unsigned int i = 1; i < gridSize; ++i){
		double testDist = DistSqd(x, y, z, grid[i]);
		if(testDist < bestDist){
			bestDist = testDist;
			bestInd = i;
		}
	}
	return bestInd;
}

//Finds the closest grid point to the grid's band maxima and returns
//its index.  In general, the smoothed energies have slightly 
//different values, so need the direction (0, 1, 2) -> (x, y, z) too.
int FindMaxima(const struct gp** grid, const int gridSize,
			   const int dir){
	int bestInd = 0;
	double bestNrg = grid[0]->nrgS[dir];
	for(unsigned int i = 1; i < gridSize; ++i){
		double testNrg = grid[i]->nrgS[dir];
		if(testNrg > bestNrg){
			bestNrg = testNrg;
			bestInd = i;
		}
	}
	return bestInd;
}

//Finds the closest grid point to the grid's band minima and returns
//its index.  In general, the smoothed energies have slightly 
//different values, so need the direction (0, 1, 2) -> (x, y, z) too.
int FindMinima(const struct gp** grid, const int gridSize,
			   const int dir){
	int bestInd = 0;
	double bestNrg = grid[0]->nrgS[dir];
	for(unsigned int i = 1; i < gridSize; ++i) {
		double testNrg = grid[i]->nrgS[dir];
		if(testNrg < bestNrg){
			bestNrg = testNrg;
			bestInd = i;
		}
	}
	return bestInd;
}

//Gives band structure along a direction (0, 1, 2 = x, y, z) 
//originating from a point in grid that has index sInd.  
//Writes results to fileName as k0 e0 eS0\nk1 e1 eS1\n ..., but not 
//necessairily in order (I move from the pivot to the low dir, then 
//the high dir)
void WriteBS(const char *fileName, const struct gp **grid, 
			 const int dir, const int sInd){
	//Init arrays
	int actSize = 0;
	double *k = malloc(MAX_SIZE*sizeof(double));
	double *e = malloc(MAX_SIZE*sizeof(double));
	double *s = malloc(MAX_SIZE*sizeof(double));

	int lo = 2*dir; int hi = 2*dir + 1;

	
	struct gp *curr = grid[sInd];
	k[0] = curr->crds[dir];
	e[0] = curr->nrg;
	s[0] = curr->nrgS[dir];
	actSize++;

	///low loop (break at low edge)
	for(; actSize < MAX_SIZE; ++actSize){
		if(curr->nbrs[lo] != NULL){
			curr = curr->nbrs[lo];
			k[actSize] = curr->crds[dir];
			e[actSize] = curr->nrg;
			s[actSize] = curr->nrgS[dir];
			continue;
		}
		break;
	}
	curr = grid[sInd];
	///high loop (break at high edge)
	for(; actSize < MAX_SIZE; ++actSize){
		if (curr->nbrs[hi] != NULL) {
			curr = curr->nbrs[hi];
			k[actSize] = curr->crds[dir];
			e[actSize] = curr->nrg;
			s[actSize] = curr->nrgS[dir];
			continue;
		}
		break;
	}

	//Write file
	FILE* outfile;
	outfile = fopen(fileName, "a");
	if (!outfile) {
		printf("WriteBS(): Unable to locate file :(\n");
		exit(1);
	}

	for(unsigned int i = 0; i < actSize; ++i){
		fprintf(outfile, "%.8f %.8f %.8f\n", k[i], e[i], s[i]);
	}

	fclose(outfile);
	free(k);
	free(e);
	free(s);
}

//Writes the energy of the entire BZ as an array of
//(kx, ky, kz, nrg, nrgSmooth[0], ... nrgSmooth[3]).  
void WriteFullBs(const char* fileName, const struct gp** grid,
				 const unsigned int gridSize){
	FILE* outfile;
	outfile = fopen(fileName, "a");
	if (!outfile) {
		printf("WriteBS(): Unable to locate file :(\n");
		exit(1);
	}
	
	for(unsigned int i = 0; i < gridSize; ++i){
		if(grid[i]->isInBz == 't'){
			fprintf(outfile, "%.8f %.8f %.8f %.8f %.8f %.8f" 
				    "%.8f\n",
					grid[i]->crds[0], grid[i]->crds[1], 
					grid[i]->crds[2], grid[i]->nrg,
					grid[i]->nrgS[0], grid[i]->nrgS[1],
					grid[i]->nrgS[2]);
		}
	}

	fclose(outfile);
}