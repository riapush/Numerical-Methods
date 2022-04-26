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

double* RK(double a, double b, double y0, double h, double(*f)(double, double)) { // with given h
	int N = (int)((b - a) / h);
	double x = a; // starting point
	double* y_vec = malloc((N + 1) * sizeof(double));
	for (int i = 0; i <= N; ++i) {
		y_vec[i] = y0;
		double k_1 = f(x, y0);
		printf("k_1 = %.5lf \n", k_1);
		double k_2 = f(x + h / 2, y0 + h * k_1 / 2.0);
		printf("k_2 = %.5lf \n", k_2);
		double k_3 = f(x + h / 2, y0 + h * k_2 / 2.0);
		printf("k_3 = %.5lf \n", k_3);
		double k_4 = f(x + h, y0 + h * k_3);
		printf("k_4 = %.5lf \n", k_4);
		y0 += (h / 6.0) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
		x += h;
	}
	return y_vec;
}

void localErr(double a, double b, double y0, double h, double(*f)(double, double)) {
	FILE* file = fopen("local.csv", "w");
	FILE* file_h = fopen("local_h.csv", "w");
	for (int i = 0; i < 15; ++i) {
		double* y = RK(a, b, y0, h, f);
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
		double* y = RK(a, b, y0, h, f);
		fprintf(file, "%.15f ", fabs(y[n] - acc(b)));
		fprintf(file_h, "%.15f ", h);
		free(y);
		h /= 2.0;
	}
	fclose(file);
	fclose(file_h);
}

void perturbation(double a, double b, double y0, double(*f)(double, double)) {
	FILE* file = fopen("pert.csv", "w");
	double h = (b - a) / 32;
	for (double i = 0.08; i < 1; i *= 2) {
		double y_per = y0 + 1.0 * i;
		double* y = RK(a, b, y_per, h, f);
		int j = 0;
		for (double x = a; x <= b; x += h) {
			fprintf(file, "%.15f ", y[j]);
			j++;
		}
		free(y);
		fprintf(file, "\n\n");
	}
	fclose(file);
}

void perturbation_small(double a, double b, double y0, double(*f)(double, double)) {
	FILE* file = fopen("pert_small.csv", "w");
	double h = (b - a) / 32;
	for (double i = 0.0001; i < 0.05; i *= 2) {
		double y_per = y0 + 1.0 * i;
		double* y = RK(a, b, y_per, h, f);
		int j = 0;
		for (double x = a; x <= b; x += h) {
			fprintf(file, "%.15f ", y[j]);
			j++;
		}
		free(y);
		fprintf(file, "\n\n");
	}
	fclose(file);
}

void perturbation_glob(double a, double b, double y0, double(*f)(double, double)) {
	FILE* glob = fopen("pert_glob.csv", "w");
	double h = (b - a) / 32;
	int N = (int)((b - a) / h);
	for (double i = 0.0001; i < 1; i *= 2) {
		double y_per = y0 + 1.0 * i;
		double* y = RK(a, b, y_per, h, f);
		int j = 0;
		fprintf(glob, "%.15f ", y[N]);
		free(y);
	}
	fclose(glob);
}

void numSolution(double a, double b, double y0, double(*f)(double, double)) {
	FILE* file = fopen("sol.csv", "w");
	double h = (b - a) / 32;
	double* y = RK(a, b, y0, h, f);
	int j = 0;
	for (double x = a; x <= b; x += h) {
		fprintf(file, "%.15f ", y[j]);
		j++;
	}
	fclose(file);
}


int main(void) {
	remove("global.csv");
	remove("global_h.csv");
	remove("pert.csv");
	remove("pert_small.csv");
	remove("pert_glob_small.csv");
	remove("pert_glob.csv");
	remove("local.csv");
	remove("local_h.csv");
	remove("sol.csv");

	double a = 0.0;
	double b = 2.0;
	double y0 = 0;
	double h = (b - a) / 2.0;
	//localErr(a, b, y0, h, dy);
	//globalErr(a, b, y0, h, dy);
	//perturbation(a, b, y0, dy);
	//perturbation_small(a, b, y0, dy);
	//perturbation_glob(a, b, y0, dy);
	//numSolution(a, b, y0, dy);

	double* y_h = RK(a, b, y0, (b - a) / 2.0, dy);
	printf("h = (b-a)/2\n y0 = %.5lf \n y1 = %.5lf \n y2 = %.5lf \n\n", y_h[0], y_h[1], y_h[2]);

	double* y_2h = RK(a, b, y0, b - a, dy);
	printf("h = (b-a)\n y0 = %.5lf \n y1 = %.5lf \n\n", y_2h[0], y_2h[1]);

	free(y_h);
	free(y_2h);
	return 0;
}