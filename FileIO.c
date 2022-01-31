
//Header file for reading, writing stuff

#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FDecs.h"



//domain: {xLo, xHi, yLo, yHi, zLo, zHi}
//gridPts: {n gridpts x, n gridpts y, n gridpts z}
//B: {{b1x, b1y, b1z}, {b2z, b2y, b2z}, {b3x, b3y, b3z}}
//smear: {smearX, smearY, smearZ}
//mus, temps, bands all follow:
// x: {x0, x1, x2, ...} of length nX
void ReadInputFile(const char *fileName,
				   double **domain, int **gridPts, 
				   double ***B, int *idwPow, 
				   double *idwRad, double **smear, 
				   unsigned int *nVTemps, double **vTemps,
				   unsigned int *nCTemps, double **cTemps, 
				   unsigned int *nVMus, double **vMus, 
				   unsigned int *nCMus, double **cMus, 
				   unsigned int *nVBands, int **vBands, 
				   unsigned int *nCBands, int **cBands){

	FILE* infile;
	infile = fopen(fileName, "r");
	if(!infile){
		printf("ReadInputFile(): Unable to locate file :(\n");
		exit(1);
	}


	char *line = malloc(sizeof(char)*LINESIZE);
	char *next = malloc(sizeof(char)*LINESIZE);
	char *trash = malloc(sizeof(char)*LINESIZE);
	char *lineCpy = line;
	char *nextCpy = next;
	char *trashCpy = trash;

	for(;;){
		fgets(line, LINESIZE, infile);

		///Termination condition
		if(strstr(line, "STOP")){
			break;
		}

		///Read cube domain
		if (strstr(line, "BEGIN DOMAIN")){
			*domain = (double*)malloc(6*sizeof(double));

			for (unsigned int i = 0; i < 3; ++i){
				fgets(line, LINESIZE, infile);
				(*domain)[2*i] = strtod(line, &next);
				(*domain)[2*i+1] = strtod(next, trash);
			}
		}

		///Read Gridpoints
		if (strstr(line, "BEGIN GRIDPOINTS")) {
			*gridPts = malloc(3*sizeof(int));

			for(unsigned int i = 0; i < 3; ++i) {
				fgets(line, LINESIZE, infile);
				(*gridPts)[i] = strtol(line, trash, BASE);
			}
		}
		
		///Read Recip Lattice
		if (strstr(line, "BEGIN LATTICE")) {
			*B = malloc(3*sizeof(double*));
			for (unsigned int i = 0; i < 3; ++i) {
				(*B)[i] = malloc(3*sizeof(double));
			}

			for (unsigned int i = 0; i < 3; ++i) {
				fgets(line, LINESIZE, infile);
				(*B)[i][0] = strtod(line, &next);
				(*B)[i][1] = strtod(next, &next);
				(*B)[i][2] = strtod(next, trash);
			}
		}
		
		///Read Inverse Distance Weighting Params
		if (strstr(line, "BEGIN IDWPARAMS")) {
			fgets(line, LINESIZE, infile);
			*idwRad = strtod(line, trash);
			fgets(line, LINESIZE, infile);
			*idwPow = strtol(line, trash, BASE);
		}

		///Read Smearing parameters
		if (strstr(line, "BEGIN SMEARPARAMS")) {
			*smear = malloc(3*sizeof(double));

			for (unsigned int i = 0; i < 3; ++i) {
				fgets(line, LINESIZE, infile);
				(*smear)[i] = strtod(line, trash);
			}
		}
		
		//Read valance temperatures
		if (strstr(line, "BEGIN VTEMPS")) {
			fgets(line, LINESIZE, infile); ///explicit or range?

			///Explicit Inputs
			if(strstr(line, "explicit")){
				fgets(line, LINESIZE, infile);
				*nVTemps = strtol(line, trash, BASE);

				*vTemps = malloc((*nVTemps)*sizeof(double));
				fgets(line, LINESIZE, infile);
				for(unsigned int i = 0; i < *nVTemps; ++i){
					(*vTemps)[i] = strtod(line, &line);
				}
			}
			///Range of inputs
			if(strstr(line, "range")){
				fgets(line, LINESIZE, infile);
				*nVTemps = strtol(line, trash, BASE);

				fgets(line, LINESIZE, infile);
				double low = strtod(line, &next);
				double high = strtod(next, trash);
				*vTemps = malloc((*nVTemps)*sizeof(double));
				FillArray_d(vTemps, *nVTemps, low, high);
			}
		}

		//Read conduction temperatures
		if (strstr(line, "BEGIN CTEMPS")) {
			fgets(line, LINESIZE, infile); ///explicit or range?

			///Explicit Inputs
			if (strstr(line, "explicit")) {
				fgets(line, LINESIZE, infile);
				*nCTemps = strtol(line, trash, BASE);

				*cTemps = malloc((*nCTemps)*sizeof(double));
				fgets(line, LINESIZE, infile);
				for (unsigned int i = 0; i < *nCTemps; ++i) {
					(*cTemps)[i] = strtod(line, &line);
				}
			}
			///Range of inputs
			if (strstr(line, "range")) {
				fgets(line, LINESIZE, infile);
				*nCTemps = strtol(line, trash, BASE);

				fgets(line, LINESIZE, infile);
				double low = strtod(line, &next);
				double high = strtod(next, trash);
				*cTemps = malloc((*nCTemps)*sizeof(double));
				FillArray_d(cTemps, *nCTemps, low, high);
			}
		}

		//Read Valance Chemical Potentials
		if (strstr(line, "BEGIN VMUS")) {
			fgets(line, LINESIZE, infile); ///explicit or range?

			///Explicit Inputs
			if (strstr(line, "explicit")) {
				fgets(line, LINESIZE, infile);
				*nVMus = strtol(line, trash, BASE);

				*vMus = malloc((*nVMus)*sizeof(double));
				fgets(line, LINESIZE, infile);
				for (unsigned int i = 0; i < *nVMus; ++i) {
					(*vMus)[i] = strtod(line, &line);
				}
			}
			///Range of inputs
			if (strstr(line, "range")) {
				fgets(line, LINESIZE, infile);
				*nVMus = strtol(line, trash, BASE);

				fgets(line, LINESIZE, infile);
				double low = strtod(line, &next);
				double high = strtod(next, trash);
				*vMus = malloc((*nVMus)*sizeof(double));
				FillArray_d(vMus, *nVMus, low, high);
			}
		}

		//Read Conduction Chemical Potentials
		if (strstr(line, "BEGIN CMUS")) {
			fgets(line, LINESIZE, infile); ///explicit or range?

			///Explicit Inputs
			if (strstr(line, "explicit")) {
				fgets(line, LINESIZE, infile);
				*nCMus = strtol(line, trash, BASE);

				*cMus = malloc((*nCMus)*sizeof(double));
				fgets(line, LINESIZE, infile);
				for (unsigned int i = 0; i < *nCMus; ++i) {
					(*cMus)[i] = strtod(line, &line);
				}
			}
			///Range of inputs
			if (strstr(line, "range")) {
				fgets(line, LINESIZE, infile);
				*nCMus = strtol(line, trash, BASE);

				fgets(line, LINESIZE, infile);
				double low = strtod(line, &next);
				double high = strtod(next, trash);
				*cMus = malloc((*nCMus)*sizeof(double));
				FillArray_d(cMus, *nCMus, low, high);
			}
		}

		//Read Valance band indices
		if (strstr(line, "BEGIN VBANDS")) {
			fgets(line, LINESIZE, infile); ///explicit or range?

			///Explicit Inputs
			if (strstr(line, "explicit")) {
				fgets(line, LINESIZE, infile);
				*nVBands = strtol(line, trash, BASE);

				*vBands = malloc((*nVBands)*sizeof(int));
				fgets(line, LINESIZE, infile);
				for (unsigned int i = 0; i < *nVBands; ++i) {
					(*vBands)[i] = strtol(line, &line, BASE);
				}
			}
			///Range of inputs
			if (strstr(line, "range")) {
				fgets(line, LINESIZE, infile);
				int low = strtol(line, &next, BASE);
				int high = strtol(next, trash, BASE);
				*nVBands = high - low + 1;
				*vBands = malloc((*nVBands)*sizeof(int));
				FillArray_i(vBands, *nVBands, low, high);
			}
		}
		
		//Read Conduction band indices
		if (strstr(line, "BEGIN CBANDS")) {
			fgets(line, LINESIZE, infile); ///explicit or range?

			///Explicit Inputs
			if (strstr(line, "explicit")) {
				fgets(line, LINESIZE, infile);
				*nCBands = strtol(line, trash, BASE);

				*cBands = malloc((*nCBands)*sizeof(int));
				fgets(line, LINESIZE, infile);
				for (unsigned int i = 0; i < *nCBands; ++i) {
					(*cBands)[i] = strtol(line, &line, BASE);
				}
			}
			///Range of inputs
			if (strstr(line, "range")) {
				fgets(line, LINESIZE, infile);
				int low = strtol(line, &next, BASE);
				int high = strtol(next, trash, BASE);
				*nCBands = high - low + 1;
				*cBands = malloc((*nCBands)*sizeof(int));
				FillArray_i(cBands, *nCBands, low, high);
			}
		}


	}

	free(lineCpy);
	free(nextCpy);
	free(trashCpy);
	fclose(infile);
	return;
}

//Reads the eigenvalues for all k points into two parallel arrays
//kps = {{k0x, k0y, k0z}, {k1x, k1y, k1z}, ...} length nKps, 3
//eigs = {{k0eig0, k0eig1, ...}, 
//        {k1eig0, k1eig1, ...}, ...} length nKps, nEigs
void ReadEigennvalueFile(const char *fileName, int *nKps,
						 double ***kps, int *nEigs, double ***eigs){
	FILE* infile;
	infile = fopen(fileName, "r");
	if (!infile) {
		printf("ReadEigenvalueFile(): Unable to locate file :(\n");
		exit(1);
	}


	char* line = malloc(sizeof(char)*LINESIZE);
	char* next = malloc(sizeof(char)*LINESIZE);
	char* trash = malloc(sizeof(char)*LINESIZE);
	char* lineCpy = line;
	char* nextCpy = next;
	char* trashCpy = trash;

	///Get number of k points and eigenvalues per k point
	fgets(line, LINESIZE, infile);
	*nKps = strtol(line, &next, BASE);
	*nEigs = strtol(next, trash, BASE);

	*kps = malloc((*nKps)*sizeof(double*));
	for (unsigned int i = 0; i < *nKps; ++i) {
		(*kps)[i] = malloc(3*sizeof(double));
	}
	*eigs = malloc((*nKps)*sizeof(double*));
	for (unsigned int i = 0; i < *nKps; ++i) {
		(*eigs)[i] = malloc((*nEigs)*sizeof(double));
	}

	///Actually read the rest of the file's values
	for(unsigned int i = 0; i < *nKps; ++i){
		fgets(line, LINESIZE, infile);

		(*kps)[i][0] = strtod(line, &next); ///a
		(*kps)[i][1] = strtod(next, &next); ///b
		(*kps)[i][2] = strtod(next, trash); ///c

		for(unsigned int j = 0; j < *nEigs; ++j){
			fgets(line, LINESIZE, infile);
			(*eigs)[i][j] = strtod(line, trash);
		}
	}

	free(lineCpy);
	free(nextCpy);
	free(trashCpy);
	fclose(infile);
}

//Prints the Mij^-1 tensors as:
//type mu t xx yy zz yz zx xy
//where type = 'p' for holes, 'n' for electrons, t = temp in kelvin,
//and xx, ... etc are the inverse effective mass components scaled by
//the resting electron mass i.e. xx = m_0 / m^*_xx, etc.
void WriteResults(const char *fileName, const struct minv ***mass,
				  const int nMus, const int nTs){
	FILE *outfile;
	outfile = fopen(fileName, "a");
	if (!outfile) {
		printf("WriteResults(): Unable to locate file :(\n");
		exit(1);
	}

	struct minv *m = malloc(sizeof(struct minv*));
	fprintf(outfile, "%-5s %-6s %-9s  %-11s %-11s %-11s %-11s %-11s"
				     " %-11s\n",
					 "type", "mu", "t", "xx", "yy", "zz", "yz", "zx",
				     "xy");
	for(unsigned int i = 0; i < nMus; ++i){
		for(unsigned int j = 0; j < nTs; ++j){
			m = mass[i][j];

			fprintf(outfile, "%-5c %06.3f %09.4f  %011.8f %011.8f"
							 " %011.8f %011.8f %011.8f %011.8f\n",
							 m->type, m->mu, m->temp, 
							 MInvConv(m->numer[0], m->denom),
							 MInvConv(m->numer[1], m->denom),
							 MInvConv(m->numer[2], m->denom),
							 MInvConv(m->numer[3], m->denom),
							 MInvConv(m->numer[4], m->denom),
							 MInvConv(m->numer[5], m->denom));
		}
	}

	free(m);
	fclose(outfile);
}

//Fills arr of size n with equally spaced values from low to high
void FillArray_d(double** arr, const int n,
	const double low, const double high) {
	double diff = (high - low) / ((double)n - 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		(*arr)[i] = low + diff * (double)i;
	}
	return;
}
//Fills arr of size n with ints from low to high (inclusive)
void FillArray_i(int** arr, const int n,
	const int low, const int high) {
	for (int i = 0; i < n; ++i) {
		(*arr)[i] = low + i;
	}
	return;
}

