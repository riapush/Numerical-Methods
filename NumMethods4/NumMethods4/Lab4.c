#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)

void print_matrix(double** A, int N, int M)
{
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			printf("%*.29lf ", 6, A[i][j]);
		}
		printf("\n");
	}
}


// tested, works
// creates n x m matrix
double** create_matrix(int c, int r) { // c is length of column, r is length of row
	double** A = (double**)malloc(c * sizeof(double*));
	if (A == NULL) {
		printf("Memory allocation error in create_matrix");
		return NULL;
	}
	for (int i = 0; i < c; i++) {
		A[i] = (double*)malloc(r * sizeof(double));
	}
	return A;
}


double* create_vector(int c) {
	return (double*)malloc(c * sizeof(double));
}


double** read_matrix(FILE* file_m, FILE* file_r, int* rang) {
	fscanf(file_r, "%i", rang); // read rang
	double** matrix = create_matrix(*rang, *rang); // for householders method extended matrix is used
	for (int i = 0; i < *rang; i++) {
		for (int j = 0; j < *rang; j++) {
			fscanf(file_m, "%lf; ", &matrix[i][j]);
		}
	}
	return matrix;
}

double** read_b(FILE* file_rp, int rang) {
	double** b = create_matrix(rang, 1);
	for (int i = 0; i < rang; i++) {
		fscanf(file_rp, "%lf; ", &b[i][0]);
	}
	return b;
}

// tested, works
// frees n x (any size) matrix/vector
void free_matrix(double** A, int n) { // there's n for sure
	for (int i = 0; i < n; i++) {
		free(A[i]);
	}
	free(A);
}

// multiplication of c x c matrix on c x r matrix
double** matrix_multiply(int c, int r, double** H, double** A) { //M is the rows of Matrix H, K is the columns of Matrix A, and N is the "count" of these parts that you have to add.
	double** A1 = create_matrix(c, r);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
			A1[i][j] = 0;
			for (int k = 0; k < c; k++) {
				A1[i][j] += H[i][k] * A[k][j];
			}
		}
	}
	return A1;
}

double* create_zero_vector(int n) {
	double* x = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		x[i] = 0;
	}
	return x;
}

// returns vector-row
double* vector_w(double* a, int r) {
	double norm = vector_norm_2(a, r);
	double beta = (a[0] <= 0) ? norm : -norm;
	double mu = 1 / (2 * beta*beta - 2 * beta*a[0]); // we check that norm != 0 in householders method

	double* w = (double*)create_vector(r);
	w[0] = mu * (a[0] - beta);
	for (int i = 0; i < r; i++) {
		w[i] = mu * a[i];
	}
	return w;
}


// householders matrix is square.
double** householder_matrix(double* a, int ñ) {

	double* w = vector_w(a, ñ);
	double** H = create_e_matrix(ñ, ñ);

	for (int i = 0; i < ñ; i++) {
		for (int j = 0; j < ñ; j++) {
			H[i][j] -= 2 * (w[i] * w[j]);
		}
	}

	free(w);
	return H;
}

double** enlarge_matrix(double** m, int n) {
	double** new_m = create_matrix(n + 1, n + 1);
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			if (i == 0 && j == 0) {
				new_m[i][j] = 1;
			}
			else if (i == 0 || j == 0) {
				new_m[i][j] = 0;
			}
			else {
				new_m[i][j] = m[i - 1][j - 1];
			}
		}
	}
	free_matrix(m, n);
	return new_m;
}


void householders_qr_decomposition(double** A, int c, double*** Q, double*** R) {
	double** A1 = create_matrix(c, c);
	for (int i = 0; i < c; i++) {
		double** A_t = transposed_matrix(A, c, c);

		double* u = (A_t[i] + i);

		double* v = create_vector(c - i);
		if (u[0] < 0) {
			v[0] = vector_norm_2(u, c - i);
		}
		else {
			v[0] = -vector_norm_2(u, c - i);
		}
		for (int j = 1; j < c - i; j++) {
			v[j] = 0;
		}

		double** w = create_matrix(1, c);
		for (int j = 0; j < c - i; j++) {
			w[0][j] = u[j] - v[j];
		}
		for (int j = c - i; j < c; j++) {
			w[0][j] = 0;
		}
		double w_norm = vector_norm_2(w[0], c - i);
		if (w_norm == 0) {
			continue;
		}
		double** P = create_e_matrix(c - i, c - i);
		for (int m = 0; m < c - i; m++) {
			for (int n = 0; n < c - i; n++) {
				P[m][n] -= (2 / (w_norm*w_norm)) * w[0][m] * w[0][n];
			}
		}
		int rang_diff = c - i;
		while (rang_diff != c) {
			P = enlarge_matrix(P, rang_diff);
			rang_diff++;
		}
		if (i = 0) {
			(*Q) = P;
		}
		else {
			double** Q1 = matrix_multiply(c, c, (*Q), P);
			free_matrix((*Q), c);
			(*Q) = Q1;
		}
		A_n(c, i, w, A);


		free_matrix(A_t, c);
		free(v);
		free(w);
	}
	(*R) = A;
	free_matrix(A, c);
}



double* qr_doubleshift_decomposition(double** A, int n, int eps) {
	double** Q = NULL;
	double** R = NULL;

	double** P = create_e_matrix(c - i, c - i);
	for (int m = 0; m < c - i; m++) {
		for (int n = 0; n < c - i; n++) {
			P[m][n] -= (2 / (w_norm*w_norm)) * w[0][m] * w[0][n];
		}
	}
	int rang_diff = c - i;
	while (rang_diff != c) {
		P = enlarge_matrix(P, rang_diff);
		rang_diff++;
	}
}

