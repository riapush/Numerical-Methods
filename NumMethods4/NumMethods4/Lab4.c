#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#pragma warning(disable : 4996)


void print_matrix(double** A, int N, int M)
{
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			printf("%*.2lf& ", 6, A[i][j]);
		}
		printf("\\\\\n");
	}
}

// creates n x m matrix
double** create_matrix(int n) { // c is length of column, r is length of row
	double** A = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; i++) {
		A[i] = (double*)malloc(n * sizeof(double));
	}
	return A;
}

double** create_w(int c, int r) {
	double** A = (double**)malloc(c * sizeof(double*));
	for (int i = 0; i < c; i++) {
		A[i] = (double*)malloc(r * sizeof(double));
	}
	return A;
}

double** read_matrix(FILE* file_m, int rang) {
	double** matrix = create_matrix(rang);
	for (int i = 0; i < rang; i++) {
		for (int j = 0; j < rang; j++) {
			fscanf(file_m, "%lf; ", &matrix[i][j]);
		}
	}
	return matrix;
}

double* create_vector(int c) {
	return (double*)malloc(c * sizeof(double));
}

void free_matrix(double** A, int n) { // there's n for sure
	for (int i = 0; i < n; i++) {
		free(A[i]);
	}
	free(A);
}

double** transposed_matrix(double** a, int n) { // transpose matrix
	double** a_t = create_matrix(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a_t[j][i] = a[i][j];
		}
	}
	return a_t;
}

double** create_E_matrix(int c) {
	double** A = create_matrix(c);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < c; j++) {
			if (i == j)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}
	return A;
}


double vector_norm_2(double* A, int n) {
	double norm = 0;
	for (int i = 0; i < n; i++) {
		norm += A[i] * A[i];
	}
	return sqrt(norm);
}



double** matrix_multiply(int n, double** H, double** A) {
	double** A1 = create_matrix(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A1[i][j] = 0;
			for (int k = 0; k < n; k++) {
				A1[i][j] += H[i][k] * A[k][j];
			}
		}
	}
	return A1;
}


double** hessenberg(double** A, int n) {
	double** B = A;
	double summa;
	double s;
	int sign;
	double m;
	for (int i = 0; i < n - 2; i++) {
		summa = 0;

		for (int j = i + 1; j < n; j++) {
			summa += B[j][i] * B[j][i];
		}
		if (-B[i + 1][i] > 0) {
			sign = 1;
		}
		else {
			sign = -1;
		}
		s = sign * sqrt(summa);

		m = 1 / sqrt(2 * s * (s - B[i + 1][i]));

		double* w = create_vector(n);
		for (int j = 0; j < n; j++) {
			if (j <= i) {
				w[j] = 0;
			}
		}

		w[i + 1] = m * (B[i + 1][i] - s);

		for (int p = i + 2; p < n; p++) {
			w[p] = m * B[p][i];
		}

		double** P = create_E_matrix(n);
		for (int v = 0; v < n; v++) {
			for (int k = 0; k < n; k++) {
				P[v][k] -= 2 * w[v] * w[k];
			}
		}
		B = matrix_multiply(n, P, B);
		B = matrix_multiply(n, B, P);

	}
	return B;
}

void givens(double** B, double*** Q, double*** R, int n) {
	double t;
	double cos;
	double sin;
	for (int j = 0; j < n - 1; j++) {
		t = B[j][j] / B[j + 1][j];
		cos = 1 / sqrt(1 + t * t);
		sin = t * cos;
		double** G = create_E_matrix(n);
		G[j][j] = sin;
		G[j + 1][j] = -cos;
		G[j + 1][j + 1] = sin;
		G[j][j + 1] = cos;
		B = matrix_multiply(n, G, B);
		if (j == 0) {
			(*Q) = transposed_matrix(G, n);
		}
		else {
			double** G_t = transposed_matrix(G, n);
			double** Q1 = matrix_multiply(n, (*Q), G_t);
			free_matrix((*Q), n);
			free_matrix(G_t, n);
			free_matrix(G, n);
			(*Q) = Q1;
		}
	}
	(*R) = B;
}


double find_max_underdiagonal_element(double** matrix, int n) {
	double max = 0;
	for (int j = 0; j < (n - 1); j++) {
		for (int i = (j + 1); i < n; i++) {
			if (fabs(matrix[i][j] > max))
				max = fabs(matrix[i][j]);
		}
	}
	return max;
}

int stop_iteration(double** matrix, int n, double eps) {
	if (find_max_underdiagonal_element(matrix, n) < eps)
		return 1;
	return 0;
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

double** enlarge_matrix(double** m, int n) {
	double** new_m = create_matrix(n + 1);
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

void householders_method(double** A, double*** Q, double*** R, int c) {
	double** A1 = create_matrix(c);
	for (int i = 0; i < c; i++) {
		double** A_t = transposed_matrix(A, c);

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

		double** w = create_w(1, c);
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
		double** P = create_E_matrix(c - i);
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
		//printf("P =\n");
		//print_matrix(P, c, c);
		A_n(c, i, w, A);

		if (i == 0) {
			(*Q) = P;
		}
		else {
			double** Q1 = matrix_multiply(c, (*Q), P);
			free_matrix((*Q), c);
			(*Q) = Q1;
		}

		free_matrix(A_t, c);
		free(v);
		free(w);
	}
	(*R) = A;
}

double* qr(double** A, int n, double eps) {
	double** Q = create_matrix(n);
	double** R = create_matrix(n);
	double** A1;
	int iter = 0;
	do {
		householders_method(A, &Q, &R, n);
		A1 = matrix_multiply(n, R, Q);
		free_matrix(A, n);
		A = A1;
		iter++;
	} while (!stop_iteration(A, n, eps));
	FILE* file_iter = fopen("iterations.csv", "a");
	fprintf(file_iter, "%i;", iter);
	fclose(file_iter);

	double* solution = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		solution[i] = A[i][i];
	}

	//printf("\n");
	return solution;
}


double* qr_hessenberg(double** A, double*** eig, int n, double eps) {
	double** Q = NULL;
	double** R = NULL;
	double** A1;
	int iter = 0;
	double** B = hessenberg(A, n);
	do {
		givens(B, &Q, &R, n);
		A1 = matrix_multiply(n, R, Q);
		B = A1;
		iter++;
		if (iter == 1) {
			(*eig) = Q;
		}
		else {
			double** eig1 = matrix_multiply(n, Q, (*eig));
			(*eig) = eig1;
		}
	} while (!stop_iteration(B, n, eps));
	FILE* file_iter = fopen("iterations_hes.csv", "a");
	fprintf(file_iter, "%i;", iter);
	fclose(file_iter);

	double* solution = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		solution[i] = B[i][i];
	}

	return solution;
}

//double* qr_hessenberg(double** A, double*** eig, int n, double eps) {
//	double** Q = NULL;
//	double** R = NULL;
//	double** A1;
//	int iter = 0;
//	double** B = NULL;
//	do {
//		B = hessenberg(A, n);
//		givens(B, &Q, &R, n);
//		A1 = matrix_multiply(n, R, Q);
//		A = A1;
//		iter++;
//		if (iter == 1) {
//			(*eig) = Q;
//		}
//		else {
//			double** eig1 = matrix_multiply(n, Q, (*eig));
//			(*eig) = eig1;
//		}
//	} while (!stop_iteration(A, n, eps));
//	FILE* file_iter = fopen("iterations_hes.csv", "a");
//	fprintf(file_iter, "%i;", iter);
//	fclose(file_iter);
//
//	double* solution = (double*)malloc(n * sizeof(double));
//	for (int i = 0; i < n; i++) {
//		solution[i] = B[i][i];
//	}
//
//	return solution;
//}

void matrix_difference(double*** matrix1, double** matrix2, int n) {
	for (int i = 0; i < n; i++) {
		(*matrix1)[i][i] -= matrix2[i][i];
	}
}

void matrix_sum(double*** matrix1, double** matrix2, int n) {
	for (int i = 0; i < n; i++) {
		(*matrix1)[i][i] += matrix2[i][i];
	}
}

double* qr_iterations_shift(double** A, int n, double eps) {
	int iter = 0;
	int n_curr = n;
	int e = 0;
	double* solution = (double*)malloc(n * sizeof(double));
	double** Q = NULL;
	double** R = NULL;
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
		double** E_a = create_E_matrix(n_curr);
		for (int i = 0; i < n_curr; i++) {
			E_a[i][i] *= A[n_curr - 1][n_curr - 1];
		}
		matrix_difference(&A, E_a, n_curr);
		double** B = hessenberg(A, n_curr);
		givens(B, &Q, &R, n_curr);
		double** A1 = matrix_multiply(n_curr, R, Q);
		matrix_sum(&A1, E_a, n_curr);
		A = A1;
		double approx_zero = eps + fabs(A[n_curr - 1][n_curr - 1]);
		double vector_norm = vector_norm_2(A[n_curr - 1], n_curr);
		if (vector_norm < approx_zero && vector_norm >= fabs(A[n_curr - 1][n_curr - 1])) {
			solution[e] = A[n_curr - 1][n_curr - 1];
			e++;
			n_curr--;
		}
		iter++;
	} while (!(stop_iteration(A, n_curr, eps) && n_curr != 0));
	FILE* file_iter = fopen("iterations_shift.csv", "a");
	fprintf(file_iter, "%i;", iter);
	fclose(file_iter);

	for (int i = 0; e + i < n; i++) {
		solution[e + i] = A[i][i];
	}
	return solution;
}

int main(void) {
	//double** m = create_matrix(4);
	//m[0][0] = 5;
	//m[1][0] = 3;
	//m[2][0] = -4;
	//m[3][0] = 1;

	//m[0][1] = 1;
	//m[1][1] = 0;
	//m[2][1] = -1;
	//m[3][1] = 1;

	//m[0][2] = -3;
	//m[1][2] = -2;
	//m[2][2] = 1;
	//m[3][2] = 1;

	//m[0][3] = 1;
	//m[1][3] = 1;
	//m[2][3] = 1;
	//m[3][3] = 1;

	//print_matrix(m, 4, 4);
	//double** eig = NULL;
	//double* solution = qr_hessenberg(m, &eig, 4, pow(10,-2));
	//printf("qr_shift = %.3lf\n %.3lf\n %.3lf\n %.3lf\n\n", solution[0], solution[1], solution[2], solution[3]);
	//for (int m = 0; m < 4; m++) {
	//	for (int n = 0; n < 4; n++) {
	//		printf("%.7lf;", eig[m][n]);
	//	}
	//	printf("\n");
	//}

	//m = create_matrix(4);
	//m[0][0] = 5;
	//m[1][0] = 3;
	//m[2][0] = -4;
	//m[3][0] = 1;

	//m[0][1] = 1;
	//m[1][1] = 0;
	//m[2][1] = -1;
	//m[3][1] = 1;

	//m[0][2] = -3;
	//m[1][2] = -2;
	//m[2][2] = 1;
	//m[3][2] = 1;

	//m[0][3] = 1;
	//m[1][3] = 1;
	//m[2][3] = 1;
	//m[3][3] = 1;

	//double* solution = qr_iterations_shift(m, 4, pow(10, -10));
	//printf("qr = %.3lf\n %.3lf\n %.3lf\n %.3lf\n", solution[0], solution[1], solution[2], solution[3]);

	//m = create_matrix(4);
	//m[0][0] = 5;
	//m[1][0] = 3;
	//m[2][0] = -4;
	//m[3][0] = 1;

	//m[0][1] = 1;
	//m[1][1] = 0;
	//m[2][1] = -1;
	//m[3][1] = 1;

	//m[0][2] = -3;
	//m[1][2] = -2;
	//m[2][2] = 1;
	//m[3][2] = 1;

	//m[0][3] = 1;
	//m[1][3] = 1;
	//m[2][3] = 1;
	//m[3][3] = 1;

	//double* solution = qr_hessenberg(m, 4, pow(10, -3));
	//printf("qr hes = %.3lf\n %.3lf\n %.3lf\n %.3lf\n", solution[0], solution[1], solution[2], solution[3]);


	remove("iterations.csv");
	remove("iterations_shift.csv");
	remove("iterations_hes.csv");
	FILE* fileA = fopen("matrix.csv", "r");
	for (int i = 0; i < 10; i++) {
		double** A = read_matrix(fileA, 10);
		double* solution = qr(A, 10, pow(10,-10));
		free(solution);
	}
	fclose(fileA);

	double** eig = 0;
	fileA = fopen("matrix.csv", "r");
	for (int i = 0; i < 10; i++) {
		double** A = read_matrix(fileA, 10);
		double* solution = qr_hessenberg(A, &eig, 10, pow(10, -10));
		free(solution);
	}
	fclose(fileA);

	fileA = fopen("matrix.csv", "r");
	for (int i = 0; i < 10; i++) {
		double** A = read_matrix(fileA, 10);
		double* solution = qr_iterations_shift(A, 10, pow(10, -10));
		free(solution);
	}
	fclose(fileA);


	fileA = fopen("perturbation.csv", "r");
	FILE* fileS = fopen("solution_delta.csv", "w");
	FILE* fileE = fopen("eigen_vectors.csv", "w");
	for (int i = 0; i < 10; i++) {
		double** A = read_matrix(fileA, 10);
		double* solution = qr_hessenberg(A, &eig, 10, pow(10, -10));
		for (int m = 0; m < 10; m++) {
			for (int n = 0; n < 10; n++) {
				fprintf(fileE, "%.15lf;", eig[m][n]);
			}
			fprintf(fileE, "\n");
		}
		fprintf(fileE, "\n");
		for (int j = 0; j < 10; j++) {
			fprintf(fileS, "%.15lf;", solution[j]);
		}
		fprintf(fileS, "\n");
		free(solution);
	}
	fclose(fileA);
	fclose(fileS);
	fclose(fileE);


}
