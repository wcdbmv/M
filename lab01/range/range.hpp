#ifndef LAB01_RANGE_RANGE_HPP_
#define LAB01_RANGE_RANGE_HPP_

#include <vector>

class Range {
public:
	using const_reference = std::vector<double>::const_reference;
	using const_iterator = std::vector<double>::const_iterator;

public:
	Range(double left, double right, double delta);

	[[nodiscard]] double left() const;
	[[nodiscard]] double right() const;
	[[nodiscard]] double delta() const;

	void set_left(double left);
	void set_right(double right);
	void set_delta(double delta);

	[[nodiscard]] const_reference operator[](size_t index) const;
	[[nodiscard]] const_iterator begin() const;
	[[nodiscard]] const_iterator end() const;
	[[nodiscard]] size_t size() const;

private:
	double left_;
	double right_;
	double delta_;
	std::vector<double> nodes_;

	[[nodiscard]] bool check_range(double left, double right) const;
	[[nodiscard]] bool check_delta(double delta) const;
	void calc_nodes();
};

#endif  // LAB01_RANGE_RANGE_HPP_
