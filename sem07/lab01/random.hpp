#pragma once

#include <random>
#include <unordered_map>

/*
 * today is 6 oct 2020 and Apple Clang does not support concepts
 */
template <typename T>
T randint(T min, T max) {
	static thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_int_distribution<T> distribution(min, max);
	return distribution(generator);
}

template <typename Container, typename Func>
Container generate_sequence(decltype(Container{}.size()) count, Func func) {
	Container sequence(count);
	std::generate(sequence.begin(), sequence.end(), func);
	return sequence;
}

template <typename T>
T limit(T value, T min, T max) {
	return value % (max - min + 1) + min;
}

template <template <typename> typename Container, typename T>
double serial_test(const Container<T>& values, T min, T max) {
	// size must be even
	if (values.size() % 2) {
		return 0;
	}

	const auto chi_squared_test = [](const std::unordered_map<size_t, size_t>& frequency, double p, int n) -> double {
		double chi_squared = 0;
		for (const auto& [index, y] : frequency) {
			chi_squared += y * y;
		}
		chi_squared = chi_squared / (n * p) - n;
		return chi_squared;
	};

	const auto d = max - min + 1;
	const auto d2 = d * d;
	const auto n = values.size() / 2;

	std::unordered_map<size_t, size_t> frequency;
	for (decltype(values.size()) i = 0; i < n; ++i) {
		const size_t index = static_cast<size_t>((values[2 * i] - min) * d + values[2 * i + 1] - min);
		++frequency[index];
	}

	const auto p = 1.0 / d2;
	return chi_squared_test(frequency, p, n);
}
