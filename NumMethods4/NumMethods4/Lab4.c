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
			printf("%*.5lf ", 6, A[i][j]);
		}
		printf("\n");
	}
}

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

double** read_matrix(FILE* file_m, int rang) {
	double** matrix = create_matrix(rang, rang);
	for (int i = 0; i < rang; i++) {
		for (int j = 0; j < rang; j++) {
			fscanf(file_m, "%lf; ", &matrix[i][j]);
		}
	}
	return matrix;
}

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
		//printf("s = %.3lf\n", s1);
		if (s1 == 0) {
			continue;
		}
		double mu = 1 / (sqrt(2 * s1*(s1 - A[i + 1][i])));
		//printf("mu = %.3lf\n", mu);
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
		//printf("\nH = \n");
		//print_matrix(P, c, c);
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
		//printf("\nA in hess = \n");
		//print_matrix(A, c, c);


		free_matrix(A_t, c);
		free(w);
	}
	(*R) = A;
	//printf("\nB = \n");
	//print_matrix(B, c, c);
	return B;
}

void householders_method(double** A, double*** Q, double***R, int c) {
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
		A_n(c, i, w, A);


		free_matrix(A_t, c);
		free_matrix(P, c);
		free(v);
		free(w);
	}

}

void givens_decomposition(double** A, double*** Q, double*** R, int c, int eps) {
	//printf("\n\ngivens decomp\n\n");
	for (int j = 0; j < c - 1; j++) {
		if (fabs(A[j + 1][j]) < pow(10,-eps)) {
			continue;
		}
		double t = A[j][j] / A[j + 1][j];
		double cos = 1 / sqrt(1 + t * t);
		double sin = t * cos;
		double** G = create_e_matrix(c, c);
		G[j][j] = sin;
		G[j + 1][j] = -cos;
		G[j + 1][j + 1] = sin;
		G[j][j + 1] = cos;
		//printf("\nG = \n");
		//print_matrix(G, c, c);
		double** A1 = matrix_multiply(c, c, G, A);
		free_matrix(A, c);
		A = A1;
		//printf("\nA in givens = \n");
		//print_matrix(A, c, c);
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
	//printf("\nR = \n");
	//print_matrix((*R), c, c);
	//printf("\nQ = \n");
	//print_matrix((*Q), c, c);
}


bool eigennum_found(double** A, int n, int eps) {
	double sum = 0;
	for (int i = 0; i < n-1; i++) {
		sum += fabs(A[i + 1][i]) * fabs(A[i + 1][i]);
	}
	if (sqrt(sum) >= pow(10,-eps)) {
		return false;
	}
	return true;
}

double** matrix_diff(double** A, double** B, int c) {
	double** C = create_matrix(c, c);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < c; j++) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
	return C;
}

double** matrix_sum(double** A, double** B, int c) {
	double** C = create_matrix(c, c);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < c; j++) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
	return C;
}

double* qr_iterations(double** A, int n, int eps) {
	int iter = 0;
	do {
		double** Q = NULL;
		double** R = NULL;
		double** B = hessenberg_form(A, &Q, &R, n);
		givens_decomposition(B, &Q, &R, n, eps);
		double** A1 = matrix_multiply(n, n, R, Q);
		A = A1;
		printf("\nA = \n");
		print_matrix(A, n, n);
		free_matrix(Q, n);
		free_matrix(R, n);
		iter++;
	} while (!eigennum_found(A, n, eps));
	FILE* file_iter = fopen("iterations.csv", "a");
	fprintf(file_iter, "%i;", iter);
	fclose(file_iter);

	//printf("\nA = \n");
	//print_matrix(A, n, n);
	double* solution = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		solution[i] = A[i][i];
	}
	for (int i = 0; i < n; i++) {
		printf("%lf ", solution[i]);
	}
	printf("\n");
	return solution;
}

int compare(double* x, double* y) {
	return *x - *y;
}

double* qr_iterations_shift(double** A, int n, int eps) {
	int iter = 0;
	int n_curr = n;
	int e = 0;
	double* solution = (double*)malloc(n * sizeof(double));
	do {
		if (n_curr == 2) {
			double D = (A[1][1] + A[0][0])*(A[1][1] + A[0][0]) - 4 * ((A[1][1] * A[0][0]) - A[1][0] * A[0][1]);
			solution[e] = (A[1][1] + A[0][0] + sqrt(D)) / 2;
			e++;
			solution[e] = (A[1][1] + A[0][0] - sqrt(D)) / 2;
			e++;
			break;
		}
		if (n_curr == 1) {
			solution[e] = A[0][0];
			e++;
			break;
		}
		double** Q = NULL;
		double** R = NULL;
		double** E_a = create_e_matrix(n_curr, n_curr);
		for (int i = 0; i < n_curr; i++) {
			E_a[i][i] *= A[n_curr-1][n_curr-1];
		}
		double** A_shifted = matrix_diff(A, E_a, n_curr);
		double** B = hessenberg_form(A_shifted, &Q, &R, n_curr);
		givens_decomposition(B, &Q, &R, n_curr, eps);
		double** A1 = matrix_multiply(n_curr, n_curr, R, Q);
		double** A2 = matrix_sum(A1, E_a, n_curr);
		free_matrix(E_a, n_curr);
		free_matrix(A1, n_curr);
		A = A2;
		free_matrix(Q, n_curr);
		free_matrix(R, n_curr);
		double approx_zero = pow(10, -eps) + fabs(A[n_curr - 1][n_curr - 1]);
		double vector_norm = vector_norm_2(A[n_curr - 1], n_curr);
		if (vector_norm < approx_zero && vector_norm >= fabs(A[n_curr - 1][n_curr - 1])) {
			solution[e] = A[n_curr - 1][n_curr - 1];
			e++;
			n_curr--;
		}
		//printf("\nA = \n");
		//print_matrix(A, n_curr, n_curr);
		iter++;
	} while (!(eigennum_found(A, n_curr, eps) && n_curr != 0));
	FILE* file_iter = fopen("iterations.csv", "a");
	fprintf(file_iter, "%i;", iter);
	fclose(file_iter);

	for (int i = 0; e + i < n; i++) {
		solution[e + i] = A[i][i];
	}
	qsort(solution, n, sizeof(double), (int(*) (const void *, const void *))compare);
	for (int i = 0; i < n; i++) {
		printf("%.15lf ", solution[i]);
	}
	printf("\n");
	return solution;
}


int main(void) {

	remove("accuracy.csv");
	for (int i = 1; i < 14; i++) {
		FILE* fileA = fopen("matrix_accuracy.csv", "r");
		double** A = read_matrix(fileA, 10);
		fclose(fileA);
		double* solution = qr_iterations_shift(A, 10, i);
		FILE* file_sol = fopen("accuracy.csv", "a");
		for (int j = 0; j < 10; j++) {
			fprintf(file_sol, "%.15lf;", solution[j]);
		}
		fprintf(file_sol, "\n");
		free(solution);
		fclose(file_sol);
	}
	/*remove("iterations.csv");
	FILE* fileA = fopen("matrix.csv", "r");
	for (int i = 0; i < 16; i++) {
		double** A = read_matrix(fileA, 5);
		double* solution = qr_iterations_shift(A, 5, 10);
		free(solution);
	}
	fclose(fileA);

	FILE* file_sol = fopen("solution.csv", "w");
	FILE* Aper = fopen("perturbation.csv", "r");
	for (int i = 0; i < 7; i++) {
		double** A = read_matrix(Aper, 5);
		double* solution = qr_iterations_shift(A, 5, 10);
		for (int j = 0; j < 5; j++) {
			fprintf(file_sol, "%.15lf;", solution[j]);
		}
		fprintf(file_sol, "%\n");
		free(solution);
	}
	fclose(file_sol);
	fclose(Aper);

	double** m = create_matrix(3, 3);
	m[0][0] = 5;
	m[1][0] = 3;
	m[2][0] = -4;
	m[0][1] = 1;
	m[1][1] = 0;
	m[2][1] = -1;
	m[0][2] = -3;
	m[1][2] = -2;
	m[2][2] = 1;

	print_matrix(m, 3, 3);
	double* solution = qr_iterations(m, 3, 3);
	printf("%.3lf %.3lf %.3lf", solution[0], solution[1], solution[2]);*/
}
