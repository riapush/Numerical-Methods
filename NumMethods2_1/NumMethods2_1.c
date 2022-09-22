#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)
#define DELTA 1000.0

typedef struct {
	double x;
	double y;
	double dy;
}point_t;

double func(double x) {
	return x * x - sin(10 * x);
}

double funcDer(double x) {
	return 2 * x - 10 * cos(10 * x);
}

double* chebyshevGrid(int numOfPoints, double a, double b) {
	FILE* cheb_points = fopen("cheb_points.csv", "a");
	double* roots = (double*)malloc(numOfPoints * sizeof(double));
	for (int i = 0; i < numOfPoints; i++) {
		roots[i] = cos(M_PI * (2 * i + 1) / (2 * (numOfPoints)));
		roots[i] = (b + a) / 2 + (b - a) / 2 * roots[i];
		fprintf(cheb_points, "%.15lf;", roots[i]);
	}
	fprintf(cheb_points, "\n");
	fclose(cheb_points);
	return roots;
}

double* uniformGrid(int numOfPoints, double a, double b) {
	FILE* uni_points = fopen("uni_points.csv", "a");
	double* roots = (double*)malloc(numOfPoints * sizeof(double));
	for (int i = 0; i < numOfPoints; i++) {
		roots[i] = a + i * ((b - a) / (double)(numOfPoints-1));
		fprintf(uni_points, "%.15lf;", roots[i]);
	}
	fprintf(uni_points, "\n");
	fclose(uni_points);
	return roots;
}


point_t* createTableFunc(double b, double a, int numOfPoints, double(*func)(double), double(*der)(double), double*(*grid)(int, double, double)) {
	point_t* points = (point_t*)malloc(numOfPoints * sizeof(point_t));
	double* roots = grid(numOfPoints, a, b);

	for (int i = 0; i < numOfPoints; ++i) {
		points[i].x = roots[i];
		points[i].y = func(points[i].x);
		points[i].dy = der(points[i].x);
	}

	free(roots);

	return points;
}

double Hermit(point_t* points, int numOfPoints, double x) {
	double sum = 0;
	double prod = 1;
	double res = 0;

	for (int j = 0; j < numOfPoints; ++j) {
		for (int k = 0; k < numOfPoints; ++k) {
			if (k == j)
				continue;
			sum += (x - points[j].x) / (points[j].x - points[k].x);
			prod *= (x - points[k].x) / (points[j].x - points[k].x) * (x - points[k].x) / (points[j].x - points[k].x);
		}
		res += ((x - points[j].x) * points[j].dy + (1 - 2 * sum) * points[j].y) * prod;
		sum = 0;
		prod = 1;
	}

	return res;
}


	int main(void){
		remove("cheb.csv");
		remove("uni.csv");
		remove("cheb_points.csv");
		remove("uni_points.csv");
		remove("x.csv");
		double a = -2;
		double b = 2;
		double step = (b - a) / DELTA;
		FILE* cheb = fopen("cheb.csv", "w");
		FILE* uni = fopen("uni.csv", "w");
		for (int numOfPoints = 3; numOfPoints < 25; numOfPoints++) {
			double x = a;
			point_t* points_cheb = createTableFunc(b, a, numOfPoints, func, funcDer, chebyshevGrid);
			point_t* points_uni = createTableFunc(b, a, numOfPoints, func, funcDer, uniformGrid);
			for (int i = 0; i < 100; i++) {
				if (numOfPoints == 4) {
					FILE* file_x = fopen("x.csv", "a");
					fprintf(file_x, "%.15lf;", x);
					fclose(file_x);
				}
				double res_cheb = Hermit(points_cheb, numOfPoints, x);
				fprintf(cheb, "%.15lf;", res_cheb);

				double res_uni = Hermit(points_uni, numOfPoints, x);
				fprintf(uni, "%.15lf;", res_uni);

				x += step;
			}
			fprintf(cheb, "\n");
			fprintf(uni, "\n");
		}
		fclose(cheb);
		fclose(uni);

		FILE* comparison = fopen("comparison_1.csv", "a");
		double x = a;
		for (int i = 0; i < 1000; i++) {
			point_t* points_uni = createTableFunc(b, a, 15, func, funcDer, uniformGrid);
			double res_cheb = Hermit(points_uni, 15, x);
			fprintf(comparison, "%.15lf;", res_cheb);

			x += step;
		}
		fclose(comparison);
	}