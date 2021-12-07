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


double** hessenberg_form(double** A, double*** Q, double*** R, int c) { // reduction to the Hessenberg form by Householders method
	double** B = NULL;
	double** A1 = create_matrix(c, c);
	for (int i = 0; i < c - 2; i++) {
		double** A_t = transposed_matrix(A, c, c);

		double* u = (A_t[i] + i + 1);
		double s1 = (A[i + 1][i] >= 0) ? -1 : 1;
		s1 *= vector_norm_2(u, c - i - 1);
		printf("s = %.3lf\n", s1);
		if (s1 == 0) {
			continue;
		}
		double mu = 1 / (sqrt(2 * s1*(s1 - A[i + 1][i])));
		printf("mu = %.3lf\n", mu);
		double** w = create_matrix(1, c);
		for (int j = 0; j < c; j++) {
			if (j <= i) {
				w[0][j] = 0;
			}
			else if (j == i + 1) {
				w[0][j] = mu * (A[j][i] - s1);
			}
			else if (j > i + 1) {
				w[0][j] = mu * A[j][i];
			}
		}

		double** P = create_e_matrix(c, c);
		for (int m = 0; m < c; m++) {
			for (int n = 0; n < c; n++) {
				P[m][n] -= 2 * w[0][m] * w[0][n];
			}
		}
		int rang_diff = c;
		while (rang_diff != c) {
			P = enlarge_matrix(P, rang_diff);
			rang_diff++;
		}
		printf("\nH = \n");
		print_matrix(P, 3, 3);
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
		if (B == NULL) {
			B = matrix_multiply(c, c, A, P);
		}
		else {
			double** B1 = matrix_multiply(c, c, P, B);
			double** B2 = matrix_multiply(c, c, B1, P);
			free_matrix(B1, c);
			free_matrix(B, c);
			B = B2;
		}
		printf("\nA = \n");
		print_matrix(A, 3, 3);


		free_matrix(A_t, c);
		//free(v);
		free(w);
	}
	(*R) = A;
	return B;
	//printf("\nR = \n");
	//print_matrix((*R), 3, 3);
	//printf("\nQ = \n");
	//print_matrix((*Q), 3, 3);
}

void givens_decomposition(double** A, double*** Q, double*** R, int c) {
	printf("\n\ngivens decomp\n\n");
	for (int j = 0; j < c - 1; j++) {
		double t = A[j][j] / A[j + 1][j];
		double cos = 1 / sqrt(1 + t * t);
		double sin = t * cos;
		double** G = create_e_matrix(c, c);
		G[j][j] = sin;
		G[j + 1][j] = -cos;
		G[j + 1][j + 1] = sin;
		G[j][j + 1] = cos;
		printf("\nG = \n");
		print_matrix(G, 3, 3);
		double** A1 = matrix_multiply(c, c, G, A);
		free_matrix(A, c);
		A = A1;
		printf("\nA = \n");
		print_matrix(A, 3, 3);
		if (j == 0) {
			(*Q) = transposed_matrix(G,c,c);
			free_matrix(G, c);
		}
		else {
			double** G_t = transposed_matrix(G, c, c);
			double** Q1 = matrix_multiply(c, c, (*Q), G_t);
			free_matrix((*Q), c);
			free_matrix(G_t, c);
			free_matrix(G, c);
			(*Q) = Q1;
		}
	}
	(*R) = A;
	printf("\nR = \n");
	print_matrix((*R), 3, 3);
	printf("\nQ = \n");
	print_matrix((*Q), 3, 3);
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

double** matrix_diff(double** A, double** B, int c) {
	double** C = create_matrix(c, c);
	for (int i = 0; i, c; i++) {
		for (int j = 0; j < c; j++) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
	return C;
}

double** matrix_sum(double** A, double** B, int c) {
	double** C = create_matrix(c, c);
	for (int i = 0; i, c; i++) {
		for (int j = 0; j < c; j++) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
	return C;
}


double* qr_iterations(double** A, int n, int eps) {
	do {
		double** Q = NULL;
		double** R = NULL;
		double** B = hessenberg_form(A, &Q, &R, n);
		givens_decomposition(B, &Q, &R, 3);
		double** A1 = matrix_multiply(n, n, R, Q);
		free_matrix(A, n);
		A = A1;
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
	//m[0][0] = m[0][2] = 1;
	//m[0][1] = -2;
	//m[1][0] = m[2][0] = 2;
	//m[1][1] = 0;
	//m[2][1] = m[2][2] = -1;
	//m[1][2] = -3;

	//m[0][0] = 5;
	//m[1][0] = -5;
	//m[2][0] = 0;
	//m[0][1] = -3;
	//m[1][1] = 2.08;
	//m[2][1] = -0.44;
	//m[0][2] = -1;
	//m[1][2] = 0.56;
	//m[2][2] = -1.08;

	//m[0][0] = 5;
	//m[1][0] = 3;
	//m[2][0] = -4;
	//m[0][1] = 1;
	//m[1][1] = 0;
	//m[2][1] = -1;
	//m[0][2] = -3;
	//m[1][2] = -2;
	//m[2][2] = 1;


	printf("m = \n");
	print_matrix(m, 3, 3);
	double** Q = NULL;
	double** R = NULL;
	//hessenberg_form(m, &Q, &R, 3);
	//givens_decomposition(m, &Q, &R, 3);
	//free_matrix(Q, 3);
	//free_matrix(R, 3);

	print_matrix(m, 3, 3);
	double* solution = qr_iterations(m, 3, 3);
	printf("\n%.3lf %.3lf %.3lf", solution[0], solution[1], solution[2]);
}
