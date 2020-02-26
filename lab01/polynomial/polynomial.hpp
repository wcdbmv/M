#ifndef LAB01_POLYNOMIAL_POLYNOMIAL_HPP_
#define LAB01_POLYNOMIAL_POLYNOMIAL_HPP_

#include <ostream>
#include <stdexcept>
#include <vector>

class Polynomial {
public:
	explicit Polynomial(const std::vector<double>& coefficients);

	[[nodiscard]] Polynomial integrate() const;
	[[nodiscard]] Polynomial integrate_from(double xi) const;

	double operator()(double x) const;

	friend Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs);
	friend Polynomial operator+(const Polynomial& polynomial, double addendum);
	friend Polynomial operator+(double addendum, const Polynomial& polynomial);
	friend Polynomial operator-(const Polynomial& lhs, const Polynomial& rhs);
	friend Polynomial operator-(const Polynomial& minuend, double subtrahend);
	friend Polynomial operator-(double minuend, const Polynomial& subtrahend);
	friend Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs);
	friend Polynomial operator*(const Polynomial& polynomial, double factor);
	friend Polynomial operator*(double factor, const Polynomial& polynomial);
	friend Polynomial operator^(const Polynomial& base, size_t power);

	friend std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial);

private:
	size_t degree_;
	std::vector<double> coefficients_;
};

#endif  // LAB01_POLYNOMIAL_POLYNOMIAL_HPP_
