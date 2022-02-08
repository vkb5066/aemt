#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include "FDecs.h"
#include "extern//kdtree-master//kdtree.h"


const char INFILE_LOC[255] = "../in.aem";
const char EIGEN_LOC[255] = "../eig.aem";
const char OUTFILE_LOC[255] = "../out.aem";

const int N_RECIP_IMGS = 3;
unsigned int N_THREADS = 1;

int main(int argc, char *argv[]){
	//Input initializations, CLIs
	const char *infile = INFILE_LOC; const char *eigfile = EIGEN_LOC;
	const char *outfile = OUTFILE_LOC;
	for(int i = 0; i < argc; i++){
		if(argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][1] == 'I'))
			infile = argv[i + 1];
		if(argv[i][0] == '-' && (argv[i][1] == 'e' || argv[i][1] == 'E'))
			eigfile = argv[i + 1];
		if(argv[i][0] == '-' && (argv[i][1] == 'o' || argv[i][1] == 'O'))
			outfile = argv[i + 1];
		if(argv[i][0] == '-' && (argv[i][1] == 'n' || argv[i][1] == 'N'))
			N_THREADS = atoi(argv[i + 1]);
	}
	double *domain; int *gridPoints; double **B;
	int idwPow; double idwRad; double *smearParams;
	unsigned int nVTemps; double *vTemps; 
	unsigned int nCTemps; double *cTemps;
	unsigned int nVMus; double *vMus; 
	unsigned int nCMus; double *cMus;
	unsigned int nVBands; int *vBands; 
	unsigned int nCBands; int *cBands;
	unsigned int nKps; double **kps; unsigned int nEigs; double **eigs;


	//Read inputs
	printf("Reading inputs ... ");
	ReadInputFile(infile, &domain, &gridPoints, &B, 
				  &idwPow, &idwRad, &smearParams, 
				  &nVTemps, &vTemps, &nCTemps, &cTemps,
				  &nVMus, &vMus, &nCMus, &cMus,
				  &nVBands, &vBands, &nCBands, &cBands);
	ReadEigennvalueFile(eigfile, &nKps, &kps, &nEigs, &eigs);
	printf("done\n");


	//Get kd tree of recip lat points and k points
	printf("Setting kd trees ... ");
	double **recipPoints = GetRecipLatPts(N_RECIP_IMGS, B);
	void *kdNodeLat = kd_create(3); 
	void *kdNodeKpt = kd_create(3);
	SetKdTree(kdNodeLat, IntPow_i(2*N_RECIP_IMGS + 1, 3), recipPoints);
	SetKdTree(kdNodeKpt, nKps, kps);
	///No longer need explicit k-points or lattice 
	FreeArr_nxm_d(&kps, nKps);
	FreeArr_nxm_d(&recipPoints, IntPow_i(2*N_RECIP_IMGS + 1, 3));
	FreeArr_nxm_d(&B, 3);
	printf("done\n");


	//Get the derivative / integration grid + connection tables
	printf("Initializing eval grid ... ");
	int **kpConTable; double **kpIDisTable; int *kcdEntryLens;
	unsigned int gridSize;
	struct gp **grid = GiveInitGrid(domain, gridPoints, kdNodeLat, 
									kdNodeKpt, idwRad, &gridSize, 
									&kpConTable, &kpIDisTable, 
									&kcdEntryLens);
	unsigned int **sizes;
	SplitArr(grid, gridSize, N_THREADS, &sizes);
	///No longer need kd nodes
	kd_free(kdNodeKpt);
	kd_free(kdNodeLat);
	printf("done\n");


	//Move over all bands, mus, and temperatures:
	//Set energies, then derivs, then preform integrals
	double *diffs;
	FillDifferentals(&diffs, domain, gridPoints);
	free(domain);
	free(gridPoints);
	///First, electrons
	printf("Electrons ...\n");
	struct minv*** eMinv;
	InitMinv(&eMinv, nCMus, cMus, nCTemps, cTemps);
	for(unsigned int i = 0; i < nCBands; ++i){
		printf(" band %i / %i", i + 1, nCBands);
		SetEnergies(&grid, 0, gridSize, idwPow, idwRad, cBands[i], 
					kcdEntryLens, kpIDisTable, eigs, kpConTable);
		SetSmoothEnergies(&grid, 0, gridSize, smearParams, diffs);
		SetDerivs(&grid, 0, gridSize);

		///Uncomment to write band interpolations
		///int pivot = FindMinima(grid, gridSize, 0);
		///WriteBS("cx", grid, 0, pivot);
		///WriteBS("cy", grid, 1, pivot);
		///WriteBS("cz", grid, 2, pivot);
		/// 
		///Uncomment to write full BZ
		///WriteFullBs("fullC", grid, gridSize);

		for(unsigned int mI = 0; mI < nCMus; ++mI){
			for(unsigned int tI = 0; tI < nCTemps; ++tI){
				double *res;
				IntegrateC(grid, 0, gridSize, &res, cTemps[tI], 
						   cMus[mI], diffs);
				FillMinv(&eMinv, res, 'n', mI, tI);
				free(res);
			}
		}
		ClearStdout();
	}
	///then, holes (different distribution function)
	printf("\nHoles ...\n");
	struct minv*** hMinv;
	InitMinv(&hMinv, nVMus, vMus, nVTemps, vTemps);
	for(unsigned int i = 0; i < nVBands; ++i){
		printf(" band %i / %i", i + 1, nVBands);
		SetEnergies(&grid, 0, gridSize, idwPow, idwRad, vBands[i], 
					kcdEntryLens, kpIDisTable, eigs, kpConTable);
		SetSmoothEnergies(&grid, 0, gridSize, smearParams, diffs);
		SetDerivs(&grid, 0, gridSize);

		///Uncomment to write band structures
		///int pivot = FindMaxima(grid, gridSize, 0);
		///WriteBS("vx", grid, 0, pivot);
		///WriteBS("vy", grid, 1, pivot);
		///WriteBS("vz", grid, 2, pivot);
		///WriteDiagBS("vxy", grid, 0, 1, pivot);
		/// 
		///Uncomment to write full BZ
		///WriteFullBs("fullV", grid, gridSize);

		for(unsigned int mI = 0; mI < nVMus; ++mI){
			for(unsigned int tI = 0; tI < nVTemps; ++tI){
				double* res;
				IntegrateV(grid, 0, gridSize, &res, vTemps[tI],
						   vMus[mI], diffs);
				FillMinv(&hMinv, res, 'p', mI, tI);
				free(res);
			}
		}
		ClearStdout();
	}
	///No longer need most of the stuff on the heap
	free(smearParams);
	free(vTemps);
	free(cTemps);
	free(vMus);
	free(cMus);
	free(vBands);
	free(cBands);
	FreeArr_nxm_d(&eigs, nEigs);
	FreeArr_nxm_i(&kpConTable, gridSize);
	FreeArr_nxm_d(&kpIDisTable, gridSize);
	for(unsigned int i = 0; i < gridSize; ++i){
		free(grid[i]);
	}
	free(grid);
	free(kcdEntryLens);
	FreeArr_nxm_i(&sizes, N_THREADS);
	free(diffs);


	//Write output, clean up
	printf("\nWriting Results ... ");
	WriteResults(outfile, eMinv, nCMus, nCTemps);
	WriteResults(outfile, hMinv, nVMus, nVTemps);
	FreeMinv(&eMinv, nCMus, nCTemps);
	FreeMinv(&hMinv, nVMus, nVTemps);

	printf("done\n");
	return 0;
}
