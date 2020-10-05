#include "cauchy.hpp"
#include <cmath>
#include "polynomial/polynomial.hpp"

double f(double x, double y) {
	return x * x + y * y;
}

std::vector<Polynomial> picard_polynomials(double xi, double eta, size_t iters_count) {
	const Polynomial x2({0, 0, 1});
	std::vector<Polynomial> ys{Polynomial({eta})};
	for (size_t i = 0; i < iters_count; ++i) {
		ys.push_back(eta + (x2 + (ys[i]^2)).integrate_from(xi));
	}
	return ys;
}

std::vector<std::vector<double>> picard_iterative_method(const Range& xs, double xi, double eta, size_t iters_count) {
	const auto ys = picard_polynomials(xi, eta, iters_count);

	std::vector<std::vector<double>> result(ys.size(), std::vector<double>(xs.size()));

	for (size_t i = ys.size() - 3; i < ys.size(); ++i) {
		for (size_t j = 0; j < xs.size(); ++j) {
			result[i][j] = ys[i](xs[j]);
		}
	}

	return result;
}

std::vector<double> euler_explicit_method(const Range& xs, double xi, double eta) {
	const double h = xs.delta();

	double yn = eta;
	double xn = xi;

	std::vector<double> result(xs.size());
	for (size_t i = 0; i < xs.size(); ++i, xn += h) {
		yn = yn + h * f(xn, yn);
		result[i] = yn;
	}

	return result;
}

double get_y(double x, double y, double h, bool& ok) {
	double root = y;
	const double c = h * x * x + y;
	const double d = 1 - 4 * h * c;
	ok = d >= 0;
	if (ok) {
		root = (1 - std::sqrt(d)) / (2 * h);
	}
	return root;
}

std::vector<double> euler_implicit_method(const Range& xs, double xi, double eta) {
	const double h = xs.delta();

	double yn = eta;
	double xn = xi;

	std::vector<double> result(xs.size());
	bool ok = true;
	size_t i = 0;
	for (; i < xs.size() && ok; ++i, xn += h) {
		yn = get_y(xn, yn, h, ok);
		result[i] = yn;
	}
	for (; i < xs.size(); ++i) {
		result[i] = yn;
	}

	return result;
}
