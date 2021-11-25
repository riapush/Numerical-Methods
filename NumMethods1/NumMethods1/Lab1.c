#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#pragma warning(disable : 4996)


void bisectionMethod(double a, double b, double(*f)(double), int p,char* filename_con, char* filename_acc) {
	double c;
	int i = 0;
	double e = pow(10,-p);
	FILE* file1 = NULL;
	if (p == 5) {
		file1 = fopen(filename_con, "w");
	}
	while (fabs(b - a) > 2 * e) {
		i++;
		c = (b + a) / 2;
		if (f(a)*f(c) > 0) {
			a = c;
		}
		else {
			b = c;
		}
		if (p == 5) {
			fprintf(file1, "%i,%.15lf\n", i, c);
		}
	}
	i++;
	c = (a+b)/2;
	if (p == 5) {
		fprintf(file1, "%i,%.15lf\n", i, c);
		fclose(file1);
	}
	FILE* file2 = fopen(filename_acc, "a");
	fprintf(file2, "%i,%.15lf,%.15lf\n", i, pow(10,-p), c);
	fclose(file2);
}


void chordMethod(double a, double b, double(*f)(double), const int p, char* filename_con, char* filename_acc) { // a = x_prev, b = x on first iteration
	double x_next = 0;
	int i = 0;
	double e = pow(10, -p);
	FILE* file1 = NULL;
	if (p == 5) {
		file1 = fopen(filename_con, "w");
	}
	do {
		i++;
		x_next = b - f(b)*(b - a) / (f(b) - f(a));
		if (f(x_next) < 0) {
			a = x_next;
		}
		else if (f(x_next) > 0) {
			b = x_next;
		}
		if (p == 5) {
			fprintf(file1, "%i,%.15lf\n", i, x_next);
		}
	} while (f(x_next+e)*f(x_next-e) > 0);
	if (p == 5) {
		fclose(file1);
	}
	FILE* file2 = fopen(filename_acc, "a");
	fprintf(file2, "%i,%.15lf,%.15lf\n", i, pow(10,-p), x_next);
	fclose(file2);
}

double polynom(double x) {
	return pow(x, 6) + pow(x, 5) - 13 * pow(x, 3) - 9 * x + 2;
}


double transcendental(double x) {
	return pow(5, x) - 6 * x - 7;
}


//double f(double x) {
//	return x * x - 16;
//}
//
//void bisectionMethodf(double a, double b, double(*f)(double), int p) {
//	double c;
//	int i = 0;
//	double e = pow(10, -p);
//	while (fabs(b - a) > 2 * e) {
//		i++;
//		c = (b + a) / 2;
//		printf("\\([%f;%f]: f(a) = %f, f(b) = %f, c = %f\\)\n", a, b, f(a), f(b), c);
//		if (f(a)*f(c) > 0) {
//			a = c;
//		}
//		else {
//			b = c;
//		}
//	}
//	i++;
//	c = (a + b) / 2;
//	printf("\\(x^* = %f\\)\n", c);
//}
//
//void chordMethodf(double a, double b, double(*f)(double), const int p) { // a = x_prev, b = x on first iteration
//	double x_next = 0;
//	int i = 0;
//	double e = pow(10, -p);
//	do {
//		i++;
//		x_next = b - f(b)*(b - a) / (f(b) - f(a));
//		printf("a = %.10f\nb = %.10f\n\n", a, b);
//		if (f(x_next) < 0) {
//			a = x_next;
//		}
//		else if (f(x_next) > 0) {
//			b = x_next;
//		}
//		else {
//			printf("you fool");
//			return;
//		}
//	} while (f(x_next + e)*f(x_next - e) > 0);
//	printf("\\(x^* = %f\\)\n", x_next);
//}



int main(void) {

	// вычисляем и записываем значение х на каждом шаге, чтобы исследовать сходимость
	bisectionMethod(2,2.5,&polynom,5, "p_bisection_converge.csv", "p_bisection_accuracy.csv");
	chordMethod(2, 2.5, &polynom, 5, "p_chord_converge.csv", "p_chord_accuracy.csv");
	chordMethod(1.5, 1.9, &transcendental, 5, "t_chord_converge.csv", "t_chord_accuracy.csv");
	bisectionMethod(1.5, 1.9, &transcendental, 5, "t_bisection_converge.csv", "t_bisection_accuracy.csv");

	for (int i = 1; i <= 15; i++) {
		bisectionMethod(1.5, 2.5, &polynom, i, "p_bisection_converge.csv", "p_bisection_accuracy.csv");
		chordMethod(1.5, 2.5, &polynom, i, "p_chord_converge.csv", "p_chord_accuracy.csv");
		chordMethod(1.5, 2, &transcendental, i, "t_chord_converge.csv", "t_chord_accuracy.csv");
		bisectionMethod(1.5, 2, &transcendental, i, "t_bisection_converge.csv", "t_bisection_accuracy.csv");
	}

	for (int i = 3; i <= 13; i += 2) {
		bisectionMethod(2, i, &polynom, 5, "dont_matter.csv", "p_bisection_root.csv");
		bisectionMethod(1.4, i, &transcendental, 5, "dont_matter.csv", "t_bisection_root.csv");

		chordMethod(2, i, &polynom, 5, "dont_matter.csv", "p_chord_root.csv");
		chordMethod(1.4, i, &transcendental, 5, "dont_matter.csv", "t_chord_root.csv");
	}
	
	return 0;
}