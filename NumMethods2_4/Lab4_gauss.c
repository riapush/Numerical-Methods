#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)

static int dots = 0;

double f(double x) { return x * x - sin(10 * x); }
double f_abs(double x) { return pow(x, 1.0/3) - sin(10*x); }
double F(double x) { return 0.1*cos(10 * x) + x * x*x / 3.0; }
double I_abs(double x) { return 0.1*cos(10.*x) + 0.75*pow(x,(4 / 3.0)); }

double x(double x_k, double a, double b) {
	return (a + b) / 2.0 + (b - a) / 2.0*x_k;
}

double quadGaussFormula(double a, double b, double(*f)(double)) {
	double A[4] = { (18 - sqrt(30)) / 36, (18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36 };
	double x_k[4] = { sqrt(3.0 / 7.0 + sqrt(24.0 / 245.0)), -sqrt(3.0 / 7.0 + sqrt(24.0 / 245.0)), sqrt(3.0 / 7.0 - sqrt(24.0 / 245.0)), -sqrt(3.0 / 7.0 - sqrt(24.0 / 245.0)) };

	double sum = 0;
	for (int i = 0; i < 4; ++i) {
		sum += A[i] * f(x(x_k[i], a, b));
		dots++;
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
	FILE* file_n = fopen("splits.csv", "a");
	FILE* err = fopen("err.csv", "a");
	FILE* file_vol = fopen("dots.csv", "a");
	fprintf(file_vol, "%i; ", dots);
	fprintf(file_n, "%i; ", n);
	fprintf(err, "%.15lf; ", fabs(F(b)-F(a)-I));
	fclose(file_n);
	fclose(err);
	fclose(file_vol);
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
	FILE* file_n = fopen("splits_non.csv", "a");
	FILE* file_err = fopen("err_non.csv", "a");
	FILE* file_vol = fopen("dots_non.csv", "a");
	fprintf(file_vol, "%i; ", dots);
	fprintf(file_n, "%i; ", n);
	fprintf(file_err, "%.15lf; ", fabs(I_abs(b)-I_abs(a)-I));
	fclose(file_n);
	fclose(file_err);
	fclose(file_vol);
	return I;
}


int main(void) {
	remove("splits.csv");
	remove("err.csv");
	remove("splits_non.csv");
	remove("err_non.csv");
	remove("dots.csv");
	for (int i = 1; i < 13; i++) {
		double I = quadGaussMethod(0, 2, f, pow(10, -i));
		double I1 = quadGaussMethod_nonsmooth(0, 2, f_abs, pow(10, -i));
		dots = 0;
	}

	//printf("%.3lf ", quadGaussMethod(0, 2, f, pow(10, -6)));

	return 0;
}