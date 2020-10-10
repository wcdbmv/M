#pragma once

#include <random>

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
