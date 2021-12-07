#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#pragma warning(disable : 4996)

void print_matrix(double** A, int N, int M)
{
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			printf("%.3lf ", A[i][j]);
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

double** transposed_matrix(double** a, int c, int r) {
	double** a_t = create_matrix(r, c);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
			a_t[j][i] = a[i][j];
		}
	}
	return a_t;
}

double** create_e_matrix(int c, int r) {
	double** A = create_matrix(c, r);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
			if (i == j)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}
	return A;
}


double vector_norm_2(double* A, int r) {
	double norm = 0;
	for (int i = 0; i < r; i++) {
		norm += A[i] * A[i];
	}
	return sqrt(norm);
}

void A_n(int n, int i0, double** w, double** A) {
	double* w_vec = *w; //one dimensional vec
	double w_norm = vector_norm_2(w_vec, n);
	double beta = 2.0 / (w_norm * w_norm);
	for (int j = 0; j + i0 < n; j++) {
		double column_base = 0;
		for (int k = 0; k + i0 < n; k++) {
			column_base += A[k + i0][j + i0] * w_vec[k];
		}
		for (int i = 0; i + i0 < n; i++) {
			A[i + i0][j + i0] = A[i + i0][j + i0] - beta * column_base * w_vec[i];
		}
	}
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


void householders_qr_decomposition(double** A, double*** Q, double*** R, int c) {
	double** A1 = create_matrix(c, c);
	double** B1 = create_matrix(c, 1);
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
		//printf("\nP = \n");
		//print_matrix(P, 3, 3);
		if (i == 0) {
			(*Q) = P;
		}
		else {
			double** Q1 = matrix_multiply(c, c, (*Q), P);
			free_matrix((*Q), c);
			(*Q) = Q1;
		}
		double** A1 = matrix_multiply(c, c, P, A);
		free_matrix(A, c);
		A = A1;
		//printf("\nA = \n");
		//print_matrix(A, 3, 3);


		free_matrix(A_t, c);
		free(v);
		free(w);
	}
	(*R) = A;
	printf("\nR = \n");
	print_matrix((*R), 3, 3);
	printf("\nQ = \n");
	print_matrix((*Q), 3, 3);
}

double* RotationsMethod(double** A, int n) {
	double* X = create_zero_vector(n);
	double cos, sin, a, b, tmp;
	for (int i = 0; i < n; i++) { //column counter
		for (int j = i + 1; j < n; j++) { //string numer counter
			a = A[i][i];
			b = A[j][i];
			cos = a / sqrt(a * a + b * b);
			sin = b / sqrt(a * a + b * b);
			// cycle for inner product k is number of a column that * on cos and sin matrix
			for (int k = i; k < n; k++) {
				tmp = A[i][k];
				A[i][k] = cos * A[i][k] + sin * A[j][k];
				A[j][k] = -sin * tmp + cos * A[j][k];
			}
		}
	}
	return X;
}


bool eigennum_found(double** A, int n, int eps) {
	double sum = 0;
	for (int i = 0; i < n-1; i++) {
		sum += A[i + 1][i] * A[i + 1][i];
	}
	if (sqrt(sum) >= pow(10,-eps)) {
		return false;
	}
	return true;
}


double* qr_iterations(double** A, int n, int eps) {
	do {
		double** Q = NULL;
		double** R = NULL;
		householders_qr_decomposition(A, &Q, &R, n);
		A = matrix_multiply(n, n, R, Q);
		printf("\nA = \n");
		print_matrix(A, 3, 3);
		free_matrix(Q, n);
		free_matrix(R, n);
	} while (!eigennum_found(A, n, eps));

	double* solution = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		solution[i] = A[i][i];
	}
	return solution;
}


int main(void) {
	double** m = create_matrix(3, 3);
	m[0][0] = 5;
	m[0][2] = -3;
	m[0][1] = 1;
	m[1][0] = 3;
	m[2][0] = -4;
	m[1][1] = 0;
	m[2][1] = -1;
	m[2][2] = 1;
	m[1][2] = -2;
	print_matrix(m, 3, 3);
	double* solution = qr_iterations(m, 3, 3);
	printf("%.3lf %.3lf %.3lf", solution[0], solution[1], solution[2]);
}
