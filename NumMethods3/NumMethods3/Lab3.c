#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
#define ESP pow(10,-3)
#pragma warning(disable:4996)


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

double* read_b(FILE* file_rp, int rang) {
	double* b = create_vector(rang);
	for (int i = 0; i < rang; i++) {
		fscanf(file_rp, "%lf; ", b + i);
	}
	return b;
}

// tested, works
// frees ñ x (any size) matrix/vector
void free_matrix(double** A, int ñ) { // there's n for sure
	for (int i = 0; i < ñ; i++) {
		free(A[i]);
	}
	free(A);
}

double* create_zero_vector(int n) {
	double* x = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		x[i] = 0;
	}
	return x;
}

double inf_norm_matrix(double** A, int n) {
	double** C = create_matrix(n, n);
	for (int i = 0; i < n; i++) {
		double curr = -A[i][i];
		for (int j = 0; j < n; j++) {
			C[i][j] = A[i][j] / curr;
		}
		C[i][i] = 0;
	}

	double max = 0;
	for (int i = 0; i < n; i++) {
		double s = 0.0;
		for (int j = 0; j < n; j++) {
			s += fabs(C[i][j]);
		}
		if (s > max)
			max = s;
	}
	free_matrix(C,n);
	return max;
}

double inf_norm_vector(double* x, int n) {
	double max = 0;
	for (int i = 0; i < n; i++){
		if (fabs(x[i]) > max)
			max = fabs(x[i]);
	}
	return max;
}


int convergence(double* xk, double* xkp, int n, double norm, double eps) {
	double* diff = create_vector(n);
	for (int i = 0; i < n; i++) {
		diff[i] = xk[i] - xkp[i];
	}
	double diff_norm = inf_norm_vector(diff, n);

	if (diff_norm < eps*(1 - norm) / norm) {
		free(diff);
		return 1;
	}
	else {
		free(diff);
		return 0;
	}
}

int diag_dominance(double** A, int n) {
	int k = 1;
	double sum;
	for (int i = 0; i < n; i++) {
		sum = 0;
		for (int j = 0; j < n; j++)
			sum += fabs(A[i][j]);
		sum -= fabs(A[i][i]);
		
		k = (sum > A[i][i]) ? 0 : 1;
	}
	return k;
}

double* gauss_seidel_iterations(int n, double** A, double* B, int eps, int* m) {
	double* x = create_zero_vector(n);
	double* p = create_vector(n);
	double norm = inf_norm_matrix(A, n);
	if (diag_dominance(A, n)) {
		do
		{
			for (int i = 0; i < n; i++)
				p[i] = x[i];
			for (int i = 0; i < n; i++)
			{
				double var = 0;
				for (int j = 0; j < n; j++)
					if (j != i)
						var += (A[i][j] * x[j]);

				x[i] = (B[i] - var) / A[i][i];
			}
			(*m)++;
			//printf("x = %.3lf %.3lf %.3lf\n", x[0], x[1], x[2]);
		} while (!convergence(x,p,n,norm,pow(10,-eps)));
	}
	return x;
}

int main(void) {
	int n = 0;
	//double** A = create_matrix(3,3);
	//A[0][0] = 4.0 / 5;
	//A[0][1] = -1.0 / 5;
	//A[0][2] = 1.0 / 5;
	//A[1][0] = -1.0 / 5;
	//A[1][1] = 3.0 / 5;
	//A[1][2] = -1.0 / 5;
	//A[2][0] = 1.0 / 5;
	//A[2][1] = -1.0 / 5;
	//A[2][2] = 5.0 / 5;
	//double* B = create_vector(3);
	//B[0] = 4.0 / 5;
	//B[1] = 1.0 / 5;
	//B[2] = 5.0 / 5;

	//double* x = gauss_seidel_iterations(3, A, B, 3, &n);

	//free(B);
	//free_matrix(A,3);

	double** A;
	double* B;
	remove("x.csv");
	remove("iter.csv");
	FILE* file_x = fopen("x.csv", "a+");
	FILE* file_y = fopen("iter.csv", "a+");

	for (int i = 1; i < 16; i++) {
		FILE* file_m = fopen("matrix.csv", "r");
		FILE* file_r = fopen("rang.csv", "r");
		FILE* file_rp = fopen("rp.csv", "r");
		A = read_matrix(file_m, file_r, &n);
		B = read_b(file_rp, n);
		int m = 0;
		double* x = gauss_seidel_iterations(n, A, B, i, &m);
		fprintf(file_y, "%.100lf;%i\n ", pow(10,-i), m);
		for (int j = 0; j < n; j++) {
			fprintf(file_x, "%.8lf;", x[j]);
		}
		fprintf(file_x, "%\n");
		m = 0;
		free(x);
		free(B);
		free_matrix(A, n, n);
		fclose(file_m);
		fclose(file_r);
		fclose(file_rp);
	}
	fclose(file_x);

	FILE* file_m = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\matrix2.csv", "r");
	FILE* file_r = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rang3.csv", "r");
	FILE* file_rp = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rp2.csv", "r");
	int m = 0;
	remove("time.csv");
	LARGE_INTEGER freq, begin, end;
	double dt;
	for (int i = 10; i <= 400; i += 10) {
		double* x;
		double** A = read_matrix(file_m, file_r, &n);
		double** B = read_b(file_rp, n);
		QueryPerformanceFrequency(&freq);
		QueryPerformanceCounter(&begin);
		for (int g = 0; g < 1000; g++) {
			x = gauss_seidel_iterations(n, A, B, 15, &m);
		}
		QueryPerformanceCounter(&end);
		dt = ((double)((end.QuadPart) - (begin.QuadPart))) / ((freq.QuadPart)*1000);
		FILE* file_t = fopen("time.csv", "a");
		fprintf(file_t, "%.20lf ", dt);
		free(x);
		fclose(file_t);
	}
	fclose(file_m);
	fclose(file_r);
	fclose(file_rp);
}