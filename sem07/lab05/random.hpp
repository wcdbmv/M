#pragma once

#include <random>

double uniform_real(double a, double b) {
	static thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_real_distribution distribution(a, b);
	return distribution(generator);
}

double normal(double mu, double sigma) {
	static thread_local std::mt19937 generator(std::random_device{}());
	std::normal_distribution distribution(mu, sigma);
	return distribution(generator);
}
