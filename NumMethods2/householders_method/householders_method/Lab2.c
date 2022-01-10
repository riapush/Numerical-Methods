#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
// frees с x (any size) matrix/vector
void free_matrix(double** A, int с) { // there's n for sure
	for (int i = 0; i < с; i++) {
		free(A[i]);
	}
	free(A);
}


// tested, works
double** transposed_matrix(double** a, int c, int r) {
	double** a_t = create_matrix(r, c);
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
			a_t[j][i] = a[i][j];
		}
	}
	return a_t;
}


// tested, works
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

double** matrix_multiply1(int c1, int r1, int r2, double** A, double** B) { //c2 = r1
	double** C = create_matrix(c1, r2);
	for (int i = 0; i < c1; i++) {
		for (int j = 0; j < r1; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < r2; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
	return C;
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

void B_n(int c, int i0, double* w, double** B) {
	double beta = 2 / vector_norm_2(w, c) / vector_norm_2(w, c);
	double dot_prod = 0;
	for (int i = 0; i + i0 < c; i++) {
		dot_prod += B[i + i0][0] * w[i];
	}
	for (int i = 0; i + i0 < c; i++) {
		B[i + i0][0] -= beta * dot_prod * w[i];
	}
}

double* reverse_gauss(double** A, double** B, int c) {
	double* x = (double*)malloc(c * sizeof(double));
	for (int i = c - 1; i >= 0; i--) {
		double sum = 0;
		if (A[i][i] == 0) {
			continue;
		}
		for (int j = i + 1; j < c; j++) {
			sum += A[i][j] * x[j];
		}
		x[i] = (B[i][0] - sum) / A[i][i];
	}
	return x;
}

// adds 1 in positions in diagonal, 0 in other positions
// only for square matrixes!
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


double* householders_method(double** A, double** B, int c) {
	double** A1 = create_matrix(c, c);
	double** B1 = create_matrix(c, 1);
	for (int i = 0; i < c; i++) {
		double** A_t = transposed_matrix(A, c, c);

		double* u = (A_t[i] + i); // iтый столбец, начиная с iтого элемента!

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
		//A1 = matrix_multiply(c, c, P, A);
		A_n(c, i, w, A);
		/*free_matrix(A, c);
		A = A1;*/
		//printf("A1 =\n");
		//for (int i = 0; i < c; i++) {
		//	for (int j = 0; j < c; j++)
		//		printf("\t%lf", A[i][j]);
		//	printf("\n");
		//}


		//B1 = matrix_multiply(c, 1, P, B);
		//free_matrix(B, c);

		B_n(c, i, w[0], B);


		free_matrix(A_t, c);
		free_matrix(P, c);
		free(v);
		free(w);
	}
	double* x = reverse_gauss(A, B, c);
	free_matrix(A, c);
	free_matrix(B, c);
	return x;
}


int main(void) {
	int n = 0;
	//double** m = create_matrix(3, 4);
	//m[0][0] = m[0][2] = 1;
	//m[0][1] = -2;
	//m[1][0] = m[2][0] = 2;
	//m[1][1] = 0;
	//m[2][1] = m[2][2] = -1;
	//m[1][2] = -3;
	//double** b = create_matrix(3, 1);
	//b[0][0] = 1;
	//b[1][0] = 8;
	//b[2][0] = 5;
	//printf("m =\n");
	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++)
	//		printf("\t%lf", m[i][j]);
	//	printf("\n");
	//}
	//double* x = householders_method(m, b, 3);
	//for (int i = 0; i < 3; i++) {
	//	printf("%lf ", x[i]);
	//}
	//free(x);


	//remove("x.csv");
	//remove("x1.csv");
	FILE* file_m = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\matrix.csv", "r");
	FILE* file_r = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rang.csv", "r");
	FILE* file_rp = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rp.csv", "r");
	//for (int i = 0; i < 10; i++) {
	//	double** A = read_matrix(file_m, file_r, &n);
	//	double** B = read_b(file_rp, n);
	//	double* x = householders_method(A, B, n);
	//	FILE* file_x = fopen("x.csv", "a");
	//	for (int j = 0; j < n; j++) {
	//		fprintf(file_x, "%.20lf ", x[j]);
	//	}
	//	fprintf(file_x, "\n");
	//	free(x);
	//	fclose(file_x);
	//}
	//fclose(file_m);
	//fclose(file_r);
	//fclose(file_rp);


	//FILE* file_r1 = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rang.csv", "r");
	//FILE* file_rp1 = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rp1.csv", "r");
	//for (int i = 0; i < 7; i++) {
	//	FILE* file_m1 = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\matrix1.csv", "r");
	//	double** A = read_matrix(file_m1, file_r1, &n);
	//	double** B = read_b(file_rp1, n);
	//	double* x = householders_method(A, B, n);
	//	FILE* file_x = fopen("x1.csv", "a");
	//	for (int j = 0; j < n; j++) {
	//		fprintf(file_x, "%.20lf ", x[j]);
	//	}
	//	fprintf(file_x, "\n");
	//	fclose(file_x);
	//	fclose(file_m1);
	//	free(x);
	//}

	//fclose(file_r1);
	//fclose(file_rp1);

	file_m = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\matrix2.csv", "r");
	file_r = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rang3.csv", "r");
	file_rp = fopen("D:\\Git\\GitHub\\Numerical-Methods\\NumMethods2\\rp2.csv", "r");

	remove("time.csv");
	for (int i = 10; i <= 400; i+=10) {
		double** A = read_matrix(file_m, file_r, &n);
		double** B = read_b(file_rp, n);
		unsigned long t_before = clock();
		double* x = householders_method(A, B, n);
		double t = ((double)clock() - (double)t_before) / CLK_TCK;
		FILE* file_t = fopen("time.csv", "a");
		fprintf(file_t, "%.20lf ", t);
		free(x);
		fclose(file_t);
	}
	fclose(file_m);
	fclose(file_r);
	fclose(file_rp);
}