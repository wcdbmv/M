#include "range.hpp"
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>

Range::Range(double left, double right, double delta) : left_(left), right_(right), delta_(delta) {
	if (!check_range(left, right)) {
		throw std::invalid_argument("Range::Range (right < left)");
	}

	if (!check_delta(delta)) {
		throw std::invalid_argument("Range::Range (delta < 0)");
	}

	calc_nodes();
}

double Range::left() const {
	return left_;
}

double Range::right() const {
	return right_;
}

double Range::delta() const {
	return delta_;
}

bool Range::check_range(double left, double right) const {
	return left <= right;
}

bool Range::check_delta(double delta) const {
	return delta > std::numeric_limits<double>::epsilon();
}

void Range::calc_nodes() {
	nodes_.resize(static_cast<size_t>(std::ceil((right_ - left_) / delta_)) + 1);

	double cur = left_;
	for (auto&& node: nodes_) {
		node = cur;
		cur += delta_;
	}
	nodes_.back() = right_;
}


void Range::set_left(double left) {
	if (!check_range(left, right_)) {
		throw std::invalid_argument("Range::set_left (right < left)");
	}

	left_ = left;

	calc_nodes();
}

void Range::set_right(double right) {
	if (!check_range(left_, right)) {
		throw std::invalid_argument("Range::set_right (right < left)");
	}

	right_ = right;

	calc_nodes();
}

void Range::set_delta(double delta) {
	if (!check_delta(delta)) {
		throw std::invalid_argument("Range::set_delta (delta < 0)");
	}

	delta_ = delta;

	calc_nodes();
}

auto Range::operator[](size_t index) const -> const_reference {
	return nodes_[index];
}

auto Range::begin() const -> const_iterator {
	return nodes_.begin();
}

auto Range::end() const -> const_iterator {
	return nodes_.end();
}

size_t Range::size() const {
	return nodes_.size();
}
