#pragma once

#include <numeric>
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
Container<double> delta_p(const Container<Container<double>>& matrix, const Container<double>& p, double dt) {
	using size_type = typename Container<double>::size_type;

	const auto n = matrix.size();
	Container<double> dp(n);
	for (size_type i = 0; i < n; ++i) {
		double curr = 0;
		for (size_type j = 0; j < n; ++j) {
			if (i == j) {
				curr += p[j] * (matrix[i][i] - std::accumulate(matrix[i].begin(), matrix[i].end(), 0.0));
			} else {
				curr += p[j] * matrix[j][i];
			}
		}
		dp[i] = curr * dt;
	}
	return dp;
}

template <template <typename> typename Container>
Container<double> solve(const Container<Container<double>>& matrix) {
	return gauss(build_coeff_matrix(matrix));
}

template <template <typename> typename Container>
Container<double> get_system_times(
		const Container<Container<double>>& matrix,
		const Container<double>& result,
		const Container<double>& p0
) {
	using size_type = typename Container<double>::size_type;
	constexpr double dt = 0.001;
	constexpr double eps = 0.01;

	const auto n = matrix.size();

	Container<double> time_result(n);
	Container<double> p = p0;
	for (auto [t, end] = std::tuple{dt, false}; !end && t < 1000; t += dt) {
		end = true;
		auto dp = delta_p(matrix, p, dt);
		for (size_type i = 0; i < n; ++i) {
			if (std::abs(result[i] - p[i]) <= eps && dp[i] <= eps) {
				if (time_result[i] == 0.0) {
					time_result[i] = t;
				}
			} else {
				end = false;
			}
			p[i] += dp[i];
		}
	}

	return time_result;
}
