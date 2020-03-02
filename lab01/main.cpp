#include <iomanip>
#include <iostream>
#include "cauchy/cauchy.hpp"

template <typename T>
T input(const std::string& prompt) {
	T value;
	std::cout << "Input " << prompt << ": ";
	std::cin >> value;
	return value;
}

int main() {
	const auto iters_count = input<size_t>("Picard's method order");
	if (iters_count <= 4) {
		std::cerr << "Error: order must be > 4\n";
		return 1;
	}
	const auto left = input<double>("min x");
	const auto right = input<double>("max x");
	const auto delta = input<double>("h");
	const Range xs(left, right, delta);

	constexpr double xi = 0.0;
	constexpr double eta = 0.0;
	const auto picard = picard_iterative_method(xs, xi, eta, iters_count);
	const auto explicit_euler = euler_explicit_method(xs, xi, eta);
	const auto implicit_euler = euler_implicit_method(xs, xi, eta);

	std::cout << "    x   │                   Picard's method                   │ Explicit method │ Implicit method\n";
	std::cout << "        │    ";
	std::cout << std::setw(2) << iters_count - 2 << " approx.   │    ";
	std::cout << std::setw(2) << iters_count - 1 << " approx.   │    ";
	std::cout << std::setw(2) << iters_count << " approx.   │                 │\n";
	std::cout << "────────┼─────────────────┼─────────────────┼─────────────────┼─────────────────┼────────────────\n";
	for (size_t i = 0; i < xs.size(); ++i) {
		std::cout << std::setw(7) << xs[i] << " │ ";
		std::cout << std::setw(15) << picard[iters_count - 2][i] << " │ ";
		std::cout << std::setw(15) << picard[iters_count - 1][i] << " │ ";
		std::cout << std::setw(15) << picard[iters_count][i] << " │ ";
		std::cout << std::setw(15) << explicit_euler[i] << " │ ";
		std::cout << std::setw(15) << implicit_euler[i] << "\n";
	}
}
