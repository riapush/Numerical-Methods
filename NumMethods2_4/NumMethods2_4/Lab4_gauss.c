#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)

double f(double x) { return x * x - sin(10 * x); }
double f_abs(double x) { return pow(x, 1.0/3) - sin(10*x); }

double x(double x_k, double a, double b) {
	return (a + b) / 2.0 + (b - a) / 2.0*x_k;
}

double quadGaussFormula(double a, double b, double(*f)(double)) {
	double A[4] = { (18 - sqrt(30)) / 36, (18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36 };
	double x_k[4] = { sqrt(3.0 / 7.0 + sqrt(24.0 / 245.0)), -sqrt(3.0 / 7.0 + sqrt(24.0 / 245.0)), sqrt(3.0 / 7.0 - sqrt(24.0 / 245.0)), -sqrt(3.0 / 7.0 - sqrt(24.0 / 245.0)) };

	double sum = 0;
	for (int i = 0; i < 4; ++i) {
		sum += A[i] * f(x(x_k[i], a, b));
	}
	return (b-a)*sum/2.0;
}

double quadGaussMethod(double a, double b, double(*f)(double), double eps) {
	int n = 1;
	double step = (b - a) / n;
	double I_prev = 0, I = 0;
	do {
		I_prev = I;
		I = 0;
		for (int i = 0; i < n; ++i) {
			I += quadGaussFormula(a + i * step, a + (i + 1)*step, f);
		}
		n *= 2;
		step = (b - a) / n;
	} while (fabs(I - I_prev) >= 255 * eps);
	FILE* file_n = fopen("iter.csv", "a");
	FILE* file_I = fopen("integral.csv", "a");
	fprintf(file_n, "%i; ", n);
	fprintf(file_I, "%.15lf; ", I);
	fclose(file_n);
	fclose(file_I);
	return I;
}

double quadGaussMethod_nonsmooth(double a, double b, double(*f)(double), double eps) {
	int n = 1;
	double step = (b - a) / n;
	double I_prev = 0, I = 0;
	do {
		I_prev = I;
		I = 0;
		for (int i = 0; i < n; ++i) {
			I += quadGaussFormula(a + i * step, a + (i + 1)*step, f);
		}
		n *= 2;
		step = (b - a) / n;
	} while (fabs(I - I_prev) >= 1.52 * eps);
	FILE* file_n = fopen("iter_non.csv", "a");
	FILE* file_I = fopen("integral_non.csv", "a");
	fprintf(file_n, "%i; ", n);
	fprintf(file_I, "%.15lf; ", I);
	fclose(file_n);
	fclose(file_I);
	return I;
}


int main(void) {
	remove("iter.csv");
	remove("integral.csv");
	remove("iter_non.csv");
	remove("integral_non.csv");
	for (int i = 1; i < 13; i++) {
		double I = quadGaussMethod(0, 2, f, pow(10, -i));
		double I1 = quadGaussMethod_nonsmooth(0, 2, f_abs, pow(10, -i));
	}

	//printf("%.3lf ", quadGaussMethod(0, 2, f, pow(10, -6)));

	return 0;
}