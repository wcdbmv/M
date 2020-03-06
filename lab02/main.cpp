#include <iostream>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <array>

// interpolation_stuff

/* constexpr std::array<std::array<double, 3>, 9> I_T0_M_table = {{
	{ 0.5,  6400, 0.40},
	{   1,  6790, 0.55},
	{   5,  7150, 1.70},
	{  10,  7270,  3.0},
	{  50,  8010, 11.0},
	{ 200,  9185, 32.0},
	{ 400, 10010, 40.0},
	{ 800, 11140, 41.0},
	{1200, 12010, 39.0}
}}; */

constexpr std::array<double, 9> I_table  = { 0.5,    1,    5,   10,   50,  200,   400,  800,   1200};
constexpr std::array<double, 9> T0_table = {6400, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010};
constexpr std::array<double, 9> M_table  = {0.40, 0.55, 1.70,  3.0, 11.0, 32.0,  40.0,  41.0,  39.0};

/* constexpr std::array<std::array<double, 2>, 11> T_sigma_table = {{
	{ 4000, 0.031},
	{ 5000,  0.27},
	{ 6000,  2.05},
	{ 7000,  6.06},
	{ 8000,  12.0},
	{ 9000,  19.9},
	{10000,  29.6},
	{11000,  41.1},
	{12000,  54.1},
	{13000,  67.7},
	{14000,  81.5}
}}; */

constexpr std::array<double, 11> T_table     = { 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000};
constexpr std::array<double, 11> sigma_table = {0.031, 0.27, 2.05, 6.06, 12.0, 19.6,  29.6,  41.1,  54.1,  67.7,  81.5};

template <typename Array>
Array logarithm(const Array& array) {
	Array result{};
	std::transform(array.begin(), array.end(), result.begin(), [](double x) { return std::log(x); });
	return result;
}

const std::array<double, 11> log_T_table = logarithm(T_table);
const std::array<double, 11> log_sigma_table = logarithm(sigma_table);

template <typename Array>
double linear_interpolation(double x, const Array& xs, const Array& ys) {
	const auto lower = std::lower_bound(xs.begin(), xs.end(), x);
	const size_t i = static_cast<size_t>(std::distance(xs.begin(), lower)) - static_cast<size_t>(lower != xs.begin());
	return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]);
}

double t0(double i) {
	return linear_interpolation(i, I_table, T0_table);
}

double m(double i) {
	return linear_interpolation(i, I_table, M_table);
}

double sigma(double t) {
	return std::exp(linear_interpolation(std::log(t), log_T_table, log_sigma_table));
}

// constants_integralCalculation

constexpr double _Tw = 2000;
constexpr double _R = 0.35;
constexpr double _Le = 12;
constexpr double _Lk = 187e-6;
constexpr double _Ck = 268e-6;
constexpr double _Rk = 0.25;
constexpr double _Uc0 = 1400;
constexpr double _I0 = 0.5;

// T(r) = T_0 + (T_w - T_0)(r / R)^m
double t(double z, double T0, double Tw, double m) {
	return T0 + (Tw - T0) * std::pow(z, m);
}

// s = h (\frac{f(x_0)}{2} + f(x_1) + \ldots + f(x_{n-1}) + \frac{f(x_n)}{2})
template <typename Function>
double trapezoidal(const Function& f, double x0, double xn, size_t n) {
	assert(x0 < xn && n > 1);

	const double h = (xn - x0) / static_cast<double>(n);
	double s = (f(x0) + f(xn)) / 2.0;
	for (size_t i = 1; i < n; ++i) {
		s += f(x0 + static_cast<double>(i)*h);
	}
	s *= h;

	return s;
}

double sigma_integral(double i) {
	return trapezoidal([i](double z) { return sigma(t(z, t0(i), _Tw, m(i))); }, 0, 1, 1024);
}

double _Rp(double i) {
	return _Le / (2 * M_PI * _R * _R * sigma_integral(i));
}

double Ucp(double I, double rp) {
	return I * rp;
}

// ord_diff_equ

double f(double I, double Uc) {
	I = std::abs(I);
	return (Uc - (_Rk + _Rp(I)) * I) / _Lk;
}

double g(double I) {
	return -std::abs(I) / _Ck;
}

std::pair<double, double> find_next_IUc(double dt, double In, double Ucn) {
	In = std::abs(In);

	const double k1 = f(In, Ucn);
	const double m1 = g(In);

	const double k2 = f(In + dt * k1 / 2.0, Ucn + dt * m1 / 2.0);
	const double m2 = g(In + dt * k1 / 2.0);

	const double k3 = f(In + dt * k2 / 2.0, Ucn + dt * m2 / 2.0);
	const double m3 = g(In + dt * k2 / 2.0);

	const double k4 = f(In + dt * k3 / 2.0, Ucn + dt * m3 / 2.0);
	const double m4 = g(In + dt * k3 / 2.0);

	const double Inp1 = In + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
	const double Ucnp1 = Ucn + dt * (m1 + 2*m2 + 2*m3 + m4) / 6.0;

	return {std::abs(Inp1), Ucnp1};
}

// main

int main() {
	// find_all_stuff
	constexpr double dt = 1e-6;
	constexpr double max_time = 800e-6;

	std::vector<double> I = {_I0};
	std::vector<double> Uc = {_Uc0};
	std::vector<double> Ucp = {_I0 * _Uc0};
	std::vector<double> Rp = {_Rp(_I0)};
	std::vector<double> time = {0};
	for (double curr_time = dt; curr_time < max_time; curr_time += dt) {
		auto [Inp1, Ucnp1] = find_next_IUc(dt, I.back(), Uc.back());
		I.push_back(Inp1);
		Uc.push_back(Ucnp1);
		Rp.push_back(_Rp(Inp1));
		Ucp.push_back(Inp1 * Rp.back());
		time.push_back(curr_time);
	}

	return 0;
}
