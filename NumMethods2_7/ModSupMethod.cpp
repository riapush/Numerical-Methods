#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <cmath>

#define SET_STREAM_PRECISION(stream) \
(stream).setf(ios::fixed);\
(stream) << setprecision(15)

using namespace std;

using vector_t = vector<double>;

struct point_t {
	double x;
	vector_t y;
	point_t(const double x_, const vector_t& y_) : x(x_), y(y_) {}
};

using grid_f = vector<point_t>;
using math_f = function<double(const double)>;
using ODE_f = function<vector_t(const double, const vector_t&)>;

struct SecOrderODE {
	math_f p, q, r, f;
	SecOrderODE(const math_f& p_, const math_f& q_, const math_f& r_, const math_f& f_) :
		p(p_), q(q_), r(r_), f(f_) {}
};

vector_t operator*(const double a, const vector_t& vec) {
	vector_t res;
	for (auto& x : vec)
		res.push_back(x * a);
	return res;
}

vector_t operator+(const vector_t& a, const vector_t& b) {
	if (a.size() != b.size())
		throw exception("Vectors adding should have same size");
	vector_t res;
	for (size_t i = 0; i < a.size(); i++)
		res.push_back(a[i] + b[i]);
	return res;
}

grid_f operator*(const double a, const grid_f& f) {
	grid_f res;
	for (auto& point : f)
		res.push_back(point_t(point.x, a * point.y));
	return res;
}

grid_f operator+(const grid_f& f, const grid_f& g) {
	grid_f res;
	for (size_t i = 0; i < f.size(); i++) {
		if (f[i].x != g[i].x)
			throw exception("GridFuncs adding should have same x grid");
		res.push_back(point_t(f[i].x, f[i].y + g[i].y));
	}
	return res;
}

ODE_f ODEFuncOf2ndOrderODE(const SecOrderODE& ode) {
	return [ode](const double x, const vector_t& y) -> vector_t {
		return { y[1], ode.f(x) / ode.p(x) - ode.r(x) / ode.p(x) * y[0] - ode.q(x) / ode.p(x) * y[1] };
	};
}

math_f zeroMathFunc() {
	return [](const double x) -> double { return 0.0; };
}

grid_f RK(const ODE_f& f, const vector_t y0, const double a, const double b, const size_t n) {
	const double h = (b - a) / (double)n;
	double x = a;
	vector_t y = y0;
	grid_f res = { point_t(a, y) };
	for (size_t i = 0; i < n; i++) {
		const vector_t k1 = f(x, y);
		const vector_t k2 = f(x + h / 2.0, y + h / 2.0 * k1);
		const vector_t k3 = f(x + h / 2.0, y + h / 2.0 * k2);
		const vector_t k4 = f(x + h, y + h * k3);
		y = y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
		x = x + h;
		res.push_back(point_t(x, y));
	}
	return res;
}

grid_f methodOfSuperposition(const SecOrderODE& ode, const double a, const double b, const vector_t& alpha, const vector_t& beta, const double A, const double B, const size_t n) {
	const vector_t u0 = { alpha[0] / (pow(alpha[0], 2) + pow(alpha[1], 2)) * A, alpha[1] / (pow(alpha[0], 2) + pow(alpha[1], 2)) * A, };
	const vector_t v0 = { alpha[1], -alpha[0]};

	grid_f u = RK(ODEFuncOf2ndOrderODE(ode), u0, a, b, n);
	grid_f v = RK(ODEFuncOf2ndOrderODE(SecOrderODE(ode.p, ode.q, ode.r, zeroMathFunc())),v0, a, b, n);
	const double c = (B - beta[0] * u[u.size() - 1].y[0] - beta[1] * u[u.size() - 1].y[1]) / (beta[0] * v[v.size() - 1].y[0] + beta[1] * v[v.size() - 1].y[1]);
	return u + c * v;
}

void writeGridFunc(const string& filename, const grid_f& gridFunc, const size_t derivativeOrder) {
	ofstream file(filename);

	SET_STREAM_PRECISION(file);
	file << "x;y" << endl;
	for (auto& point : gridFunc)
		file << point.x << ";" << point.y[derivativeOrder] << endl;
	file.close();
}

//------------------------research-------------------//

void numSolution(const string& filenameGeneric, const SecOrderODE& ode, const double a, const double b, const vector_t& alpha, const vector_t& beta, const double A, const double B, const vector<size_t>& ns, const size_t derivativeOrder) {
	for (size_t i = 0; i < ns.size(); i++) {
		printf("n = %i\n", ns[i]);
		writeGridFunc(filenameGeneric + to_string(i + 1) + ".csv", methodOfSuperposition(ode, a, b, alpha, beta, A, B, ns[i]), derivativeOrder);
	}
}

double maxError(const math_f& f, const grid_f& gridFunc, const size_t derivativeOrder) {
	double res = 0;
	for (auto& p : gridFunc) {
		double err = abs(f(p.x) - p.y[derivativeOrder]);
		if (err > res) res = err;
	}
	return res;
}

void perturbationA(const string& filename, const math_f& sol, const SecOrderODE& ode, const double a, const double b, const vector_t& alpha,
	const vector_t& beta,const double A, const double B, const size_t n, const double dA0, const size_t dASteps, const size_t derivativeOrder) {
	ofstream file(filename);
	SET_STREAM_PRECISION(file);
	file << "dy;err;" << endl;
	double dA = dA0;
	for (size_t i = 0; i < dASteps; i++) {
		file << dA << ";" << maxError( sol, methodOfSuperposition(ode, a, b, alpha, beta, A + dA, B, n), derivativeOrder) << endl;
		dA /= 10;
	}
	file.close();
}

void perturbationB(const string& filename, const math_f& sol, const SecOrderODE& ode, const double a, const double b, const vector_t& alpha,
	const vector_t& beta, const double A, const double B, const size_t n, const double dB0, const size_t dBSteps, const size_t derivativeOrder) {
	ofstream file(filename);
	SET_STREAM_PRECISION(file);
	file << "dy;err;" << endl;
	double dB = dB0;
	for (size_t i = 0; i < dBSteps; i++) {
		file << dB << ";" <<
			maxError(sol, methodOfSuperposition(ode, a, b, alpha, beta, A, B + dB, n), derivativeOrder) << endl;
		dB /= 10;
	}
	file.close();
}

void errorFromH(const string& filename, const math_f& sol, const SecOrderODE& ode, const double a, const double b, const vector_t& alpha,
	const vector_t& beta, const double A, const double B, const size_t n0, const size_t nSteps, const size_t derivativeOrder) {
	ofstream file(filename);
	SET_STREAM_PRECISION(file);
	file << "h;err;" << endl;
	size_t n = n0;
	for (size_t i = 0; i < nSteps; i++) {
		file << (b - a) / (double)n << ";" <<
			maxError(sol, methodOfSuperposition(ode, a, b, alpha, beta, A, B, n), derivativeOrder) << endl;
		n *= 2;
	}
	file.close();
}



int main() {
	const math_f acc_solution = [](const double x) -> double { return exp(x); };
	const math_f acc_solution_der = [](const double x) -> double { return exp(x); }; // 4th variant 
	const math_f p = [](const double x) -> double { return 1.0; };
	const math_f q = [](const double x) -> double { return 1.0 + pow(sin(x), 2); };
	const math_f r = [](const double x) -> double { return pow(cos(x), 2); };
	const math_f f = [](const double x) -> double { return 3.0 * exp(x); };
	const double a = 0.0;
	const double b = 1.0;
	const vector_t alpha = { 1.0, 1.0 };
	const vector_t beta = { 1.0, 0.0 };
	const double A = alpha[0] * acc_solution(a) + alpha[1] * acc_solution_der(a);
	const double B = beta[0] * acc_solution(b) + beta[1] * acc_solution_der(b);
	const size_t N = 4;
	const SecOrderODE ode(p, q, r, f);
	const vector<size_t> ns = { 4, 10 };
	const size_t n0 = 1;
	const size_t nSteps = 10;
	const size_t n = 10;
	const double dA0 = A;
	const size_t dASteps = 10;
	const double dB0 = B;
	const size_t dBSteps = 10;

	//--------------------------------------------------------

	numSolution("solution", ode, a, b, alpha, beta, A, B, ns, 0);
	//errorFromH("error_h.csv", acc_solution, ode, a, b, alpha, beta, A, B, n0, nSteps, 0);
	//perturbationA("per_A.csv", acc_solution, ode, a, b, alpha, beta, A, B, n, dA0, dASteps, 0);
	//perturbationB("per_B.csv", acc_solution, ode, a, b, alpha, beta, A, B, n, dB0, dBSteps, 0);

	return 0;
}
