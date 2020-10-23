#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

template <template <typename> typename Container>
bool is_square_matrix(const Container<Container<double>>& matrix) {
	for (const auto& row : matrix) {
		if (!row.size() == matrix.size()) {
			return false;
		}
	}
	return true;
}

template <template <typename> typename Container>
decltype(auto) max_abs_element_in_column(Container<Container<double>>& matrix, typename Container<double>::size_type col) {
	return std::max_element(matrix.begin(), matrix.end(), [col](const Container<double>& lhs, const Container<double>& rhs) {
		return fabs(lhs[col]) < fabs(rhs[col]);
	});
}

template <template <typename> typename Container>
void make_zeros(Container<Container<double>>& matrix, typename Container<double>::size_type col) {
	using size_type = typename Container<double>::size_type;

	for (size_type i = col + 1; i < matrix.size(); ++i) {
		const auto d = -matrix[i][col] / matrix[col][col];
		matrix[i][col] = 0;
		for (size_type j = col + 1; j < matrix[0].size(); ++j) {
			matrix[i][j] += d * matrix[col][j];
		}
	}
}

template <template <typename> typename Container>
Container<double> gauss(const Container<Container<double>>& A, const Container<double>&B) {
	const auto n = A.size();

	if (!n || A[0].size() != n || B.size() != n) {
		throw std::invalid_argument("Inappropriate dimensions");
	}

	using size_type = typename Container<double>::size_type;

	auto G = A;
	for (size_type k = 0; k < n; ++k) {
		auto max_abs_row = max_abs_element_in_column(G, k);
		std::swap(*max_abs_row, G[k]);
		make_zeros(G, k);

		if (fabs(G[k][k]) < std::numeric_limits<double>::epsilon()) {
			throw std::invalid_argument("hmm");
		}
	}

	Container<double> X(n);
	for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
		const auto j = static_cast<size_type>(i);
		X[j] = G[j][n] / G[j][j];
		for (int k = i - 1; k >= 0; --k) {
			const auto q = static_cast<size_type>(k);
			G[q][n] -= G[q][j] * X[j];
		}
	}

	return X;
}
