#include <iomanip>
#include <iostream>
#include "cauchy/cauchy.hpp"

int main() {
	const Range xs(0, 2.5, 0.1);
	const auto jis = picard_iterative_method(xs, 0, 0, 7);
	for (size_t i = 0; i < jis[0].size(); ++i) {
		std::cout << std::setw(6) << xs[i];
		for (size_t j = 0; j < jis.size(); ++j) {
			std::cout  << std::setw(15) << jis[j][i];
		}
		std::cout << '\n';
	}
}
