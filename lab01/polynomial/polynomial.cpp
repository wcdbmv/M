#include "polynomial.hpp"
#include <experimental/iterator>

Polynomial::Polynomial(const std::vector<double>& coefficients) :
	degree_(coefficients.size() - 1),
	coefficients_(coefficients) {
	if (coefficients.empty()) {
		throw std::invalid_argument("Coefficients must contain at least 1 element");
	}
}

Polynomial Polynomial::integrate() const {
	std::vector<double> result_coefficients(coefficients_.size() + 1);
	for (size_t i = 1; i < result_coefficients.size(); ++i) {
		result_coefficients[i] = coefficients_[i - 1] / static_cast<double>(i);
	}
	return Polynomial(result_coefficients);
}

double Polynomial::operator()(double x) const {
	double result = coefficients_[degree_];
	for (size_t i = degree_; i--;) {
		result = result * x + coefficients_[i];
	}
	return result;
}

Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs) {
	const auto plus = [](const Polynomial& min, const Polynomial& max) -> Polynomial {
		Polynomial result(max);
		for (size_t i = 0; i < min.coefficients_.size(); ++i) {
			result.coefficients_[i] += min.coefficients_[i];
		}
		return result;
	};

	return lhs.degree_ < rhs.degree_ ? plus(lhs, rhs) : plus(rhs, lhs);
}

Polynomial operator+(const Polynomial& polynomial, double addendum) {
	Polynomial result(polynomial);
	result.coefficients_[0] += addendum;
	return result;
}

Polynomial operator+(double addendum, const Polynomial& polynomial) {
	return polynomial + addendum;
}

Polynomial operator-(const Polynomial& minuend, const Polynomial& subtrahend) {
	std::vector<double> result_coefficients(minuend.coefficients_);
	result_coefficients.resize(std::max(minuend.degree_, subtrahend.degree_) + 1);
	for (size_t i = 0; i <= subtrahend.degree_; ++i) {
		result_coefficients[i] -= subtrahend.coefficients_[i];
	}

	return Polynomial(result_coefficients);
}

Polynomial operator-(const Polynomial& minuend, double subtrahend) {
	Polynomial result(minuend);
	result.coefficients_[0] -= subtrahend;
	return result;
}

Polynomial operator-(double minuend, const Polynomial& subtrahend) {
	Polynomial result(subtrahend);
	for (auto& coefficient: result.coefficients_) {
		coefficient = -coefficient;
	}
	result.coefficients_[0] += minuend;
	return result;
}

Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs) {
	std::vector<double> result_coefficients(lhs.degree_ + rhs.degree_ + 1);
	for (size_t i = 0; i <= lhs.degree_; ++i) {
		for (size_t j = 0; j <= rhs.degree_; ++j) {
			result_coefficients[i + j] += lhs.coefficients_[i] * rhs.coefficients_[j];
		}
	}
	return Polynomial(result_coefficients);
}

Polynomial operator*(const Polynomial& polynomial, double factor) noexcept {
	Polynomial result(polynomial);
	for (auto& coefficient: result.coefficients_) {
		coefficient *= factor;
	}
	return result;
}

Polynomial operator*(double factor, const Polynomial& polynomial) noexcept {
	return polynomial * factor;
}

template <typename Container>
void print_container(std::ostream& os, const Container& container) {
	std::copy(container.begin(), container.end(), std::experimental::make_ostream_joiner(os, " "));
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vector) {
	print_container(os, vector);
	return os;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial) {
	return os << polynomial.coefficients_;
}
