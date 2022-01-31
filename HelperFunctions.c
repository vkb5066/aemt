
//Header file for doing misc messy processes
#include <stdio.h>
#include "FDecs.h"

//Set KD tree, including the index of each  point as extra data to be
//accessed later
void SetKdTree(void **kdNode, const int arrSize, const double **arr){
	for(unsigned int i = 0; i < arrSize; ++i){
		kd_insert3(kdNode, arr[i][0], arr[i][1], arr[i][2], i);
	}
}

//Free n x m int array
void FreeArr_nxm_i(int ***arr, const int n){
	for(unsigned i = 0; i < n; ++i) {
		free((*arr)[i]);
	}
	free(*arr);
}

//Free n x m double array
void FreeArr_nxm_d(double ***arr, const int n){
	for(unsigned i = 0; i < n; ++i){
		free((*arr)[i]);
	}
	free(*arr);
}

//Fill array of size n with n doubles between lo and hi (inclusive)
void Linspace(double **arr, const int n, const double lo, 
		      const double hi){
	double diff = (hi - lo)/((double)n - 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		(*arr)[i] = lo + diff * (double)i;
	}
	return;
}

//Split array arr (size arrLen) into n (roughly) equally sized 
//arrays, stored in splits with sizes stored in sizes
void SplitArr(const struct gp **arr, const int arrLen,
			  const int n, unsigned int ***sizes){
	*sizes = (int**)malloc(n*sizeof(int*));

	int chunkSize = arrLen / n;
	int extra = arrLen - chunkSize*n;
	for(int i = 0, beg = 0, end = chunkSize; beg < arrLen; 
		beg = end, end = beg + chunkSize, ++i){
		if(extra){
			end++;
			extra--;
		}
		(*sizes)[i] = malloc(2*sizeof(int));
		(*sizes)[i][0] = beg;
		(*sizes)[i][1] = end;
	}
}



//Creates an 2 arrays of m inverse entries, accessed uniquely as
//arr[muI][tempI] = mInv entry of size n mus x n temps
//muI and tempI are the indices of the mus and temperatures supplied
//to the function
void InitMinv(struct minv ****arr, const int nMus, 
			   const double *mus, const int nTemps, 
			   const double *temps){
	*arr = malloc(nMus*sizeof(struct minv**));

	for(unsigned int i = 0; i < nMus; ++i){
		(*arr)[i] = malloc(nTemps*sizeof(struct minv*));

		for(unsigned int j = 0; j < nTemps; ++j){
			struct minv *p = malloc(sizeof(struct minv));
			for(unsigned int k = 0; k < 6; k++){
				p->numer[k] = 0.0;
			}
			p->denom = 0.0;

			p->mu = mus[i];
			p->temp = temps[j];

			(*arr)[i][j] = p;
		}
	}
}

void FillMinv(struct minv ****arr, const double *res, const char type,
			  const int iMu, const int iT){
	(*arr)[iMu][iT]->type = type;
	for(unsigned int i = 0; i < 6; ++i){
		(*arr)[iMu][iT]->numer[i] += res[i];
	}
	(*arr)[iMu][iT]->denom += res[6];
}

void FreeMinv(struct minv ****arr, const int iMu, const int iT){
	for(unsigned int i = 0; i < iMu; ++i){
//		for(unsigned int j = 0; j < iT; ++j){ ///this loop seems to 
//			free((*arr)[i][j]);               ///break things - not
//		}                                     ///sure why
		free((*arr)[i]);
	}
	free(*arr);
}

void ClearStdout(){
	for(unsigned int i = 0; i < LINESIZE; ++i){
		printf("\b");
	}
}