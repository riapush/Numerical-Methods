#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)

static int dots = 0;
static int dots_opt = 0;

double f(double x) { return x * x - sin(10 * x); }
double f_abs(double x) { return pow(x, 1.0 / 3) - sin(10 * x); }
double F(double x) { return 0.1*cos(10 * x) + x * x*x / 3.0; }
double I_abs(double x) { return 0.1*cos(10.*x) + 0.75*pow(x, (4 / 3.0)); }

double trapezoidalFormula(double a, double b, unsigned long long int n, double(*f)(double)) {
	double sum = f(a) + f(b);
	dots += 2;
	double h = (b - a) / (double)n;
	for (int i = 1; i < n; ++i) {
		double x = a + i * h;
		sum += 2 * f(x);
		dots++;
	}
	return 0.5*h*sum;
}

double trapezoidalFormulaOpt(double a, double b, unsigned long long int n, double(*f)(double), double I) {
	if (n == 1) {
		dots_opt += 2;
		return 0.5*(b - a)*(f(b) + f(a));
	}
	double added_sum = 0;
	double x = a + ((b - a) / (double)n); // we skip a bc it's already counted in I, same for b.
	for (int i = 0; i < n / 2; ++i) {
		added_sum += f(x) * (b - a) / (double)n;
		x += 2 * (b - a) / (double)n;
		dots_opt++;
	}
	return I * 0.5 + added_sum;
}

double trapezoidalMethod(double a, double b, double eps, double(*f)(double), char* filename) {
	int n = 1;
	double prev = 0;
	double curr = trapezoidalFormula(a, b, n, f);
	do {
		n *= 2;
		prev = curr;
		curr = trapezoidalFormula(a, b, n, f);
	} while (fabs(curr - prev) > 3 * eps);
	FILE* iter = fopen(filename, "a");
	fprintf(iter, "%i; ", dots);
	fclose(iter);
	FILE* splits = fopen("splits.csv", "a");
	fprintf(splits, "%i; ", n);
	fclose(splits);

	return curr;
}

double trapezoidalMethodNon(double a, double b, double eps, double(*f)(double), char* filename) {
	int n = 1;
	double prev = 0;
	double curr = trapezoidalFormula(a, b, n, f);
	do {
		n *= 2;
		prev = curr;
		curr = trapezoidalFormula(a, b, n, f);
	} while (fabs(curr - prev) > eps);
	FILE* iter = fopen(filename, "a");
	fprintf(iter, "%i; ", dots);
	fclose(iter);
	FILE* splits = fopen("splits_non.csv", "a");
	fprintf(splits, "%i; ", n);
	fclose(splits);

	return curr;
}

double trapezoidalMethodOpt(double a, double b, double eps, double(*f)(double), char* filename) {
	unsigned long long int n = 1;
	double prev = 0, curr;
	do {
		if (n != 1) prev = curr;
		curr = trapezoidalFormulaOpt(a, b, n, f, prev);
		n *= 2;
	} while (fabs(curr - prev) > 3 * eps);
	FILE* iter = fopen(filename, "a");
	fprintf(iter, "%i; ", dots_opt);
	fclose(iter);

	return curr;
}

int main(void) {
	remove("err.csv");
	remove("err_non.csv");
	remove("iter.csv");
	remove("iter_opt.csv");
	remove("splits.csv");
	remove("splits_non.csv");
	remove("dots_non.csv");
	remove("dots.csv");
	remove("dots_opt.csv");
	remove("dots_opt_non.csv");
	FILE* err = fopen("err.csv", "w");
	FILE* err_non = fopen("err_non.csv", "w");
	for (int n = 1; n < 13; ++n) {
		fprintf(err_non, "%.15lf; ", fabs(I_abs(2) - I_abs(0) - trapezoidalMethodNon(0, 2, pow(10, -n), f_abs, "dots_non.csv")));
		fprintf(err, "%.15lf; ", fabs(F(2)-F(0)-trapezoidalMethod(0, 2, pow(10, -n), f, "dots.csv")));

		dots = 0;
		dots_opt = 0;
	}
	fclose(err);
	fclose(err_non);

	//FILE* err_check = fopen("err_check.csv", "w"); //checking order of non-smooth func
	//FILE* h = fopen("h.csv", "w");
	//for (int n = 1; n < 1025; n *= 2) {
	//	fprintf(err_check, "%.15lf; ", fabs(I_abs(2) - I_abs(0) - trapezoidalFormula(0, 2, n, f_abs)));
	//	fprintf(h, "%.15lf; ", 2.0/n);
	//}
}