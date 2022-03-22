#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)

static dots = 0;
static dots_opt = 0;

double f(double x) {
	return x * x - sin(10 * x);
}

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

double trapezoidalMethod(double a, double b, double eps, double(*f)(double)) {
	unsigned long long int n = 1;
	double prev = 0;
	double curr = trapezoidalFormula(a, b, n, f);
	do {
		n *= 2;
		prev = curr;
		curr = trapezoidalFormula(a, b, n, f);
	} while (fabs(curr - prev) > 3 * eps);
	FILE* iter = fopen("iter.csv", "a");
	fprintf(iter, "%i; ", dots);
	fclose(iter);

	return curr;
}

double trapezoidalMethodOpt(double a, double b, double eps, double(*f)(double)) {
	unsigned long long int n = 1;
	double prev = 0, curr;
	do {
		if (n != 1) prev = curr;
		curr = trapezoidalFormulaOpt(a, b, n, f, prev);
		n *= 2;
	} while (fabs(curr - prev) > 3 * eps);
	FILE* iter = fopen("iter_opt.csv", "a");
	fprintf(iter, "%i; ", dots_opt);
	fclose(iter);

	return curr;
}

int main(void) {
	remove("err.csv");
	remove("err_opt.csv");
	remove("iter.csv");
	remove("iter_opt.csv");
	FILE* err = fopen("err.csv", "w");
	FILE* err_opt = fopen("err_opt.csv", "w");
	for (int n = 1; n < 14; ++n) {
		fprintf(err, "%.15lf; ", trapezoidalMethod(0, 2, pow(10, -n), f));
		fprintf(err_opt, "%.15lf; ", trapezoidalMethodOpt(0, 2, pow(10, -n), f));
		dots = 0;
		dots_opt = 0;
	}
	fclose(err);
	fclose(err_opt);
}