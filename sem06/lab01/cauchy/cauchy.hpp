#ifndef LAB01_CAUCHY_CAUCHY_HPP_
#define LAB01_CAUCHY_CAUCHY_HPP_

#include "range/range.hpp"

// u' = f(x, u)
// u^{(0)}(\xi) = \eta

// u' = x^2 + y^2
// u^{(0)}(0) = 0

std::vector<std::vector<double>> picard_iterative_method(const Range& xs, double xi, double eta, size_t iters_count);
std::vector<double> euler_explicit_method(const Range& xs, double xi, double eta);
std::vector<double> euler_implicit_method(const Range& xs, double xi, double eta);

#endif  // LAB01_CAUCHY_CAUCHY_HPP_
