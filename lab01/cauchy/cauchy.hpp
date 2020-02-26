#ifndef LAB01_CAUCHY_CAUCHY_HPP_
#define LAB01_CAUCHY_CAUCHY_HPP_

#include "range/range.hpp"

// y' = f(x, y)
// y^{(0)}(\xi) = \eta

// y' = x^2 + y^2
// y^{(0)}(0) = 0

std::vector<std::vector<double>> picard_iterative_method(const Range& xs, double xi, double eta, size_t iters_count);
std::vector<double> euler_explicit_method(const Range& xs, double xi, double eta);
std::vector<double> euler_implicit_method(const Range& xs, double xi, double eta);

#endif  // LAB01_CAUCHY_CAUCHY_HPP_
