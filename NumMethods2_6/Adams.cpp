#include <stdio.h>
#include <math.h>
#include <iostream> 
#pragma warning(disable : 4996)

using namespace std;

//double dy(double x, double y) {
//	return (pow(x*x + 1, 2) + 2 * x*y) / (x*x + 1);;
//}
//
//double acc(double x) { // accurate solution 
//	return x * (x*x + 1);
//}

static int volume = 0;

double dy(double x, double y) {
	volume++;
	return ((4 * x + 2 * y) / (2 * x + 1));
}

double acc(double x) {
	return ((2 * x + 1)*log(2 * x + 1) + 1);
}

double* Adams(double a, double b, double y0, double h , double(*dy)(double, double)) {
	int n = (int)(b - a) / h;
	double* ans = (double*)malloc(sizeof(double)*n);
	ans[0] = y0; // initial value
	double x = a;

	for (int i = 1; i <= 3; i++) { // RK 4th order
		//coefs
		//double k1 = h * dy(x, ans[i - 1]);
		//double k2 = h * dy(x + h / 2.0, ans[i - 1] + k1 / 2.0);
		//double k3 = h * dy(x + h / 2.0, ans[i - 1] + k2 / 2.0);
		//double k4 = h * dy(x + h, ans[i - 1] + k3);

		//ans[i] = ans[i - 1] + 1 / 6.0*(k1 + 2.0*k2 + 2.0*k3 + k4); //compute next node

		x = a + i * h;
		ans[i] = acc(x);
	}

	// adams 4th order
	for (int i = 4; i <= n; i++) {
		double k1 = 55.0*dy(x, ans[i - 1]) - 59.0*dy(x - h, ans[i - 2]) + 37.0*dy(x - 2.0*h, ans[i - 3]) - 9.0*dy(x - 3.0*h, ans[i - 4]);

		ans[i] = ans[i - 1] + h / 24.0 * k1;

		x = a + i * h;
	}

	return ans;
}

void numSolution(double a, double b, double y0, double(*f)(double, double)) {
	FILE* file = fopen("sol.csv", "w");
	double h =((b - a) / 32.0);
	double* y = Adams(a, b, y0, h, f);
	int j = 0;
	for (double x = a; x <= b; x += h) {
		fprintf(file, "%.15f ", y[j]);
		j++;
	}
	fclose(file);
}

void localErr(double a, double b, double y0, double h, double(*f)(double, double)) {
	double* y = NULL;
	FILE* file = fopen("local.csv", "w");
	FILE* file_h = fopen("local_h.csv", "w");
	for (int i = 0; i < 12; ++i) {
		y = Adams(a, b, y0, h, f);
		fprintf(file, "%.15f ", fabs(y[4] - acc(a + 4*h)));
		fprintf(file_h, "%.15f ", h);
		//free(y);
		h /= 2.0;
	}
	fclose(file);
	fclose(file_h);
}

void globalErr(double a, double b, double y0, double h, double(*f)(double, double)) {
	FILE* file = fopen("global.csv", "w");
	FILE* file_h = fopen("global_h.csv", "w");
	for (int i = 0; i < 12; ++i) {
		int n = (int)((b - a) / h);
		double* y = Adams(a, b, y0, h, f);
		fprintf(file, "%.15f ", fabs(y[n] - acc(b)));
		fprintf(file_h, "%.15f ", h);
		//free(y);
		h /= 2.0;
	}
	fclose(file);
	fclose(file_h);
}

void perturbation1(double a, double b, double y0, double h, double(*f)(double, double)){
	FILE* file = fopen("per.csv", "a");
	FILE* file_i = fopen("per_i.csv", "w");
	int n = (b - a) / h;
	for (int i = 1; i < 12; i++) {
		double y1 = y0 + pow(10, -i);
		double* y = Adams(a, b, y1, h, f);
		for (int j = 0; j <= n; j++) {
			fprintf(file, "%.15f ", y[j]);
		}
		fprintf(file, "\n");
		fprintf(file_i, "%.15f ", pow(10,-i));
	}
	fprintf(file, "\n");
	fclose(file);
	fclose(file_i);
}

void perturbation2(double a, double b, double y0, double h, double(*f)(double, double)) {
	FILE* file = fopen("per2.csv", "a");
	FILE* file_i = fopen("per_i.csv", "w");
	int n = (b - a) / h;
	for (int i = 1; i < 12; i++) {
		double y1 = y0 + pow(10, -i);
		double* y = Adams(a, b, y1, h, f);
		for (int j = 0; j <= n; j++) {
			fprintf(file, "%.15f ", y[j]);
		}
		fprintf(file, "\n");
		fprintf(file_i, "%.15f ", pow(10, -i));
	}
	fprintf(file, "\n");
	fclose(file);
	fclose(file_i);
}

void perturbation3(double a, double b, double y0, double h, double(*f)(double, double)) {
	FILE* file = fopen("per3.csv", "a");
	FILE* file_i = fopen("per_i.csv", "w");
	int n = (b - a) / h;
	for (int i = 1; i < 12; i++) {
		double y1 = y0 + pow(10, -i);
		double* y = Adams(a, b, y1, h, f);
		for (int j = 0; j <= n; j++) {
			fprintf(file, "%.15f ", y[j]);
		}
		fprintf(file, "\n");
		fprintf(file_i, "%.15f ", pow(10, -i));
	}
	fprintf(file, "\n");
	fclose(file);
	fclose(file_i);
}

int main() {
	remove("global.csv");
	remove("global_h.csv");
	remove("local.csv");
	remove("local_h.csv");
	remove("sol.csv");
	remove("per.csv");
	remove("per_i.csv");
	//double a = 0.0;
	//double b = 2.0;
	double a = 0.0;
	double b = 4;
	int n = 8;
	double h = (b - a) / n;
	//double y0 = 0.0;
	double y0 = 1.0;


	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(15);

	localErr(a, b, y0, h, dy);
	globalErr(a, b, y0, h, dy);
	numSolution(a, b, y0, dy);
	perturbation1(a, b, y0, 0.5, dy);
	perturbation2(a, b, y0, 0.25, dy);
	perturbation3(a, b, y0, 0.125, dy);

	FILE* file = fopen("vol.csv", "w");
	for (int i = 8; i <= 128; i *= 2) {
		volume = 0;
		double h1 = (b - a) / i;
		double* y = Adams(a, b, y0, h1, dy);
		fprintf(file, "%i ", volume);
	}
	fclose(file);



	return 0;
}
