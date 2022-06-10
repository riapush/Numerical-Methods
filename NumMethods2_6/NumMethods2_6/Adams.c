#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#pragma warning(disable : 4996)

double dy(double x, double y) {
	return (pow(x*x + 1, 2) + 2 * x*y) / (x*x + 1); // 15 variant
}

double acc(double x) { // accurate solution 
	return x * (x*x + 1);
}

double* Adams(double a, double b, double y0, double h, double(*f)(double, double)) {
	// we find first 4 values with RK and then use Adams

	int N = (int)((b - a) / h);
	double x = a; // starting point
	double* y_vec = malloc((N + 1) * sizeof(double));
	for (int i = 0; i <= 3; ++i) {
		y_vec[i] = y0;
		double k_1 = f(x, y0);
		double k_2 = f(x + h / 2, y0 + h * k_1 / 2.0);
		double k_3 = f(x + h / 2, y0 + h * k_2 / 2.0);
		double k_4 = f(x + h, y0 + h * k_3);
		y0 += (h / 6.0) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
		x += h;
	}
	y_vec[3] = y0;

	for (int i = 4; i <= N; ++i) {
		y_vec[i] = y_vec[i - 1] + h/24.0 * (55.0*f(x - h, y_vec[i - 1]) - 
			59.0*f(x - 2 * h, y_vec[i - 2]) + 37.0*f(x - 3 * h, y_vec[i - 3]) - 
			9.0*f(x - 4 * h, y_vec[i - 4]));

		x += h;
	}

	return y_vec;
}

void numSolution(double a, double b, double y0, double(*f)(double, double)) {
	FILE* file = fopen("sol.csv", "w");
	double h = (b - a) / 32;
	double* y = Adams(a, b, y0, h, f);
	int j = 0;
	for (double x = a; x <= b; x += h) {
		fprintf(file, "%.15f ", y[j]);
		j++;
	}
	fclose(file);
}

void localErr(double a, double b, double y0, double h, double(*f)(double, double)) {
	FILE* file = fopen("local.csv", "w");
	FILE* file_h = fopen("local_h.csv", "w");
	for (int i = 0; i < 15; ++i) {
		double* y = Adams(a, b, y0, h, f);
		fprintf(file, "%.15f ", fabs(y[1] - acc(a + h)));
		fprintf(file_h, "%.15f ", h);
		free(y);
		h /= 2.0;
	}
	fclose(file);
	fclose(file_h);
}

void globalErr(double a, double b, double y0, double h, double(*f)(double, double)) {
	FILE* file = fopen("global.csv", "w");
	FILE* file_h = fopen("global_h.csv", "w");
	for (int i = 0; i < 15; ++i) {
		int n = (int)((b - a) / h);
		double* y = Adams(a, b, y0, h, f);
		fprintf(file, "%.15f ", fabs(y[n] - acc(b)));
		fprintf(file_h, "%.15f ", h);
		free(y);
		h /= 2.0;
	}
	fclose(file);
	fclose(file_h);
}

int main(void) {
	remove("global.csv");
	remove("global_h.csv");
	remove("local.csv");
	remove("local_h.csv");
	remove("sol.csv");

	double a = 0.0;
	double b = 2.0;
	double y0 = 0;
	double h = (b - a) / 4.0;
	localErr(a, b, y0, h, dy);
	globalErr(a, b, y0, h, dy);
	numSolution(a, b, y0, dy);
	return 0;
}