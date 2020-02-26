#include "cauchy.hpp"
#include "polynomial/polynomial.hpp"

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

	for (size_t i = 0; i < ys.size(); ++i) {
		for (size_t j = 0; j < xs.size(); ++j) {
			result[i][j] = ys[i](xs[j]);
		}
	}

	return result;
}
