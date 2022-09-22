#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)

typedef struct point_t {
	double x;
	double y;
}point_t;

double f(double x) {
	return x * x - sin(10 * x);
}

double* uniformGrid(int numOfPoints, double a, double b) {
	FILE* uni_points = fopen("uni_points.csv", "a");
	double* roots = (double*)malloc(numOfPoints * sizeof(double));
	for (int i = 0; i < numOfPoints; i++) {
		roots[i] = a + i * ((b - a) / (double)(numOfPoints - 1));
		fprintf(uni_points, "%.15lf;", roots[i]);
	}
	fprintf(uni_points, "\n");
	fclose(uni_points);
	return roots;
}

point_t* createTableFunc(double b, double a, int numOfPoints, double(*func)(double), double*(*grid)(int, double, double)) {
	point_t* points = (point_t*)malloc(numOfPoints * sizeof(point_t));
	double* roots = grid(numOfPoints, a, b);

	for (int i = 0; i < numOfPoints; ++i) {
		points[i].x = roots[i];
		points[i].y = func(points[i].x);
	}

	free(roots);

	return points;
}

void gaussMethod(int m, int n, double** a, double* x) {
	int i, j, k;
	for (i = 0; i < m - 1; i++) {
		for (k = i + 1; k < m; k++) {
			if (fabs(a[i][i]) < fabs(a[k][i])) { //if diagonal element(absolute value) is smaller than any elements below it
				for (j = 0; j < n; j++) { //swap the rows
					double temp = a[i][j];
					a[i][j] = a[k][j];
					a[k][j] = temp;
				}
			}
		}
		for (k = i + 1; k < m; k++) { // gauss
			double  term = a[k][i] / a[i][i];
			for (j = 0; j < n; j++) {
				a[k][j] = a[k][j] - term * a[i][j];
			}
		}

	}
	// back-substitution
	for (i = m - 1; i >= 0; i--) {
		x[i] = a[i][n - 1];
		for (j = i + 1; j < n - 1; j++) {
			x[i] = x[i] - a[i][j] * x[j];
		}
		x[i] = x[i] / a[i][i];
	}

}

double** createMatrix(int c, int r) { // c is length of column, r is length of row
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

void freeMatrix(double** A, int n) {
	for (int i = 0; i < n; i++) {
		free(A[i]);
	}
	free(A);
}

double* leastSquaresMethod(double a, double b, int degree, int num_of_points, const char* filename) {

	point_t* points = createTableFunc(b, a, num_of_points, f, uniformGrid); // array to store the x and y-axis data-points

	double* X = (double*)malloc((2 * degree + 1) * sizeof(double)); // an array of size 2 * n + 1 for storing independent components of matrix
	for (int i = 0; i < 2 * degree + 1; i++) {
		X[i] = 0;
		for (int j = 0; j < num_of_points; j++) {
			X[i] = X[i] + pow(points[j].x, i);
		}
	}
	double** B = createMatrix(degree + 1, degree + 2); 	// augmented matrix
	double* Y = (double*)malloc((degree + 1) * sizeof(double));

	for (int i = 0; i <= degree; i++) {
		Y[i] = 0;
		for (int j = 0; j < num_of_points; j++) {
			Y[i] = Y[i] + pow(points[j].x, i) * points[j].y;
		}
	}
	for (int i = 0; i <= degree; i++) {
		for (int j = 0; j <= degree; j++) {
			B[i][j] = X[i + j];
		}
	}
	
	FILE* file_stream = fopen(filename, "a");
	for (int i = 0; i <= degree; i++) {
		for (int j = 0; j <= degree; j++) {
			fprintf(file_stream, "%.15lf; ", B[i][j]);
		}
		fprintf(file_stream, "\n");
	}
	fprintf(file_stream, "\n\n");
	fclose(file_stream);

	for (int i = 0; i <= degree; i++) {
		B[i][degree + 1] = Y[i];
	}
	double* A = (double*)malloc((degree + 1) * sizeof(double));
	gaussMethod(degree + 1, degree + 2, B, A);

	free(X);
	free(Y);
	freeMatrix(B, degree + 1);
	return A;
}

int main(void) {
	remove("fixed_degree.csv");
	remove("fixed_n.csv");
	remove("uni_points.csv");
	remove("for_cond_val.csv");
	int num_of_points; //no. of data-points
	int degree; //degree of polynomial
	double a = -2;
	double b = 2;


	FILE* fixed_degree = fopen("fixed_degree.csv", "w");
	degree = 5;
	//double* A = leastSquaresMethod(a, b, degree, 100);
	//for (int i = 0; i < degree + 1; i++) {
	//printf("%lf ", A[i]);
	//}
	for (int num_of_points = 10; num_of_points <= 100; num_of_points += 1) {
		double* A = leastSquaresMethod(a, b, degree, num_of_points, "dontmatter.csv");
		for (int i = 0; i < degree + 1; ++i) {
			fprintf(fixed_degree, "%.15lf ", A[i]);
		}
		fprintf(fixed_degree, "\n");
		free(A);
	}
	fclose(fixed_degree);


	FILE* fixed_n = fopen("fixed_n.csv", "w");
	num_of_points = 50;
	for (int degree = 1; degree < 31; degree += 1) {
		double* A = leastSquaresMethod(a, b, degree, num_of_points, "for_cond_val.csv");
		for (int i = 0; i < degree + 1; ++i) {
			fprintf(fixed_n, "%.15lf; ", A[i]);
		}
		fprintf(fixed_n, "\n");
		free(A);
	}
	fclose(fixed_n);

	FILE* comparison = fopen("comparison_2.csv", "w");
	double* A = leastSquaresMethod(a, b, 14, 15, "dontmatter.csv");
	for (int i = 0; i < 15; ++i) {
		fprintf(comparison, "%.15lf; ", A[i]);
	}
	fprintf(comparison, "\n");
	free(A);
	fclose(comparison);
}