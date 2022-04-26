#include <stdio.h>

double dy(double x, double y) {
	return (pow(x*x + 1, 2) + 2 * x*y) / (x*x + 1); // 15 variant
}

double acc(double x) { // accurate solution 
	return x * (x*x + 1);
}

double* Adams(double a, double b, double y0, double h, double(*f)(double, double)) {
	// we find first 3 values with RK and then use Adams

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
		y_vec[i] = y_vec[i - 1] + h * (55.0 / 24.0*f(x - h, y_vec[i - 1]) - 59.0 / 24.0*f(x - 2 * h, y_vec[i - 2]) +
			37.0 / 24.0*f(x - 3 * h, y_vec[i - 3]) - 3.0 / 8.0*f(x - 4 * h, y_vec[i - 4]));

		x += h;
	}

	return y_vec;
}