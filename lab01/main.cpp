#include <iostream>
#include "polynomial/polynomial.hpp"

int main() {
	Polynomial p({1, 2, 3});
	std::cout << p << '\n' << p * p << '\n' << p.integrate() << '\n' << p(2.5) << '\n';
}
