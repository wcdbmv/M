#pragma once

#include <random>

template <typename T, template <typename> typename Distribution>
T randcore(T min, T max) {
	static thread_local std::mt19937 generator(std::random_device{}());
	Distribution<T> distribution(min, max);
	return distribution(generator);
}

template <typename T>
T randint(T min, T max) {
	return randcore<T, std::uniform_int_distribution>(min, max);
}

template <typename T>
T randreal(T min, T max) {
	return randcore<T, std::uniform_real_distribution>(min, max);
}
