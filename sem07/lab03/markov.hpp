#pragma once

#include "gauss.hpp"

template <template <typename> typename Container>
Container<Container<double>> build_coeff_matrix(const Container<Container<double>>& matrix) {
	using size_type = typename Container<double>::size_type;

	const auto n = matrix.size();
	Container<Container<double>> result(n, Container<double>(n + 1));
	for (size_type state = 0; state < n - 1; ++state) {
		for (size_type col = 0; col < n; ++col) {
			result[state][state] -= matrix[state][col];
		}
		for (size_type row = 0; row < n; ++row) {
			result[state][row] += matrix[row][state];
		}
	}
	for (size_type state = 0; state < n; ++state) {
		result[n - 1][state] = 1;
	}

	result[n - 1][n] = 1;

	return result;
}

template <template <typename> typename Container>
Container<double> get_system_times(const Container<Container<double>>& matrix) {
	using size_type = typename Container<double>::size_type;

	const auto n = matrix.size();
	Container<double> result = gauss(build_coeff_matrix(matrix));
	Container<double> time_result;
	for (size_type i = 0; i < n; ++i) {
		double sum1 = 0;
		double sum2 = 0;
		for (size_type j = 0; j < n; ++j) {
			sum1 += matrix[i][j];
			sum2 += matrix[j][i];
		}
		time_result.push_back((sum1 - sum2) / result[i]);
	}
	return time_result;
}
