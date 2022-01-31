

//Header file for functions involved with manuipulating the 
//reciprocal lattice

#include <stdlib.h>

//Fill an array with reciprocal lattice points
double** GetRecipLatPts(const int nPtsPerDirection,
						const double **B) {
	double b1[3] = { B[0][0], B[0][1], B[0][2] };
	double b2[3] = { B[1][0], B[1][1], B[1][2] };
	double b3[3] = { B[2][0], B[2][1], B[2][2] };

	int lineSize = (2 * nPtsPerDirection + 1);
	int planeSize = lineSize * lineSize;
	int gridSize = planeSize * lineSize;

	//Generate a line of points along b1
	double** line = malloc(lineSize*sizeof(double*));
	for (int i = 0, m = -nPtsPerDirection; m < nPtsPerDirection + 1;
		++i, ++m) {
		line[i] = malloc(3 * sizeof(double));
		for (int j = 0; j < 3; j++) {
			line[i][j] = m * b1[j];
		}
	}

	//Generate a plane of points along b2
	double** plane = malloc(planeSize * sizeof(double*));
	for (int i = 0, m = -nPtsPerDirection; m < nPtsPerDirection + 1;
		++i, ++m) {
		for (int j = 0; j < lineSize; ++j) {
			plane[i + lineSize * j] = malloc(3*sizeof(double));

			for (int k = 0; k < 3; k++) {
				plane[i + lineSize * j][k] = line[j][k] + m * b2[k];
			}
		}
	}

	//Generate a grid of points along b3
	double** grid = malloc(gridSize * sizeof(double*));
	for (int i = 0, m = -nPtsPerDirection; m < nPtsPerDirection + 1;
		++i, ++m) {
		for (int j = 0; j < planeSize; ++j) {
			grid[i + lineSize * j] = malloc(3*sizeof(double));

			for (int k = 0; k < 3; k++) {
				grid[i + lineSize * j][k] = plane[j][k] + m * b3[k];
			}
		}
	}

	//Free memory - probably unnecessary, but if a future employer is
	//looking at this, they'll want to see it :)
	for (unsigned int i = 0; i < planeSize; ++i) {
		if (i < lineSize) {
			free(line[i]);
		}
		free(plane[i]);
	}
	free(line);
	free(plane);

	return grid;
}

//Return the index of gamma from an array of 3d points
int GetGammaIndex(const void* kdNode) {
	double dum[3] = {0., 0., 0.}; ///needed for unix systems >:(
	void* res = kd_nearest3(kdNode, 0.0, 0.0, 0.0);
	int gIndex = (int)kd_res_item3(res, dum[0], dum[1], dum[2]);
	kd_res_free(res);

	return gIndex;
}