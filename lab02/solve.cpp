#include "solve.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <array>

class Solver {
public:
	explicit Solver(const Parameters& parameters) :
		R_  (parameters.R),
		Le_ (parameters.Le),
		Lk_ (parameters.Lk),
		Ck_ (parameters.Ck),
		Rk_ (parameters.Rk),
		Uc0_(parameters.Uc0),
		I0_ (parameters.I0) {}

	Dependency solve() const {
		Dependency dependency = {
			{0},
			{I0_},
			{Uc0_},
			{Rp(I0_)},
			{Ucp(I0_, Rp(I0_))},
		};

		constexpr double dt = 1e-6;
		constexpr double max_t = 800e-6;

		for (double t = dt; t < max_t; t += dt) {
			auto [Inp1, Ucnp1] = find_next_IUc(dt, dependency.I.back(), dependency.Uc.back());
			dependency.t.push_back(t);
			dependency.I.push_back(Inp1);
			dependency.Uc.push_back(Ucnp1);
			dependency.Rp.push_back(Rp(Inp1));
			dependency.Ucp.push_back(Ucp(Inp1, dependency.Rp.back()));
		}

		return dependency;
	}

private:
	const double R_;
	const double Le_;
	const double Lk_;
	const double Ck_;
	const double Rk_;
	const double Uc0_;
	const double I0_;
	static constexpr double Tw_ = 2000;

	/*
	 * Tables
	 */
	static constexpr std::array<double, 9> I_TABLE_  = { 0.5,    1,    5,   10,   50,  200,   400,   800,  1200};
	static constexpr std::array<double, 9> T0_TABLE_ = {6400, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010};
	static constexpr std::array<double, 9> m_TABLE_  = {0.40, 0.55, 1.70,  3.0, 11.0, 32.0,  40.0,  41.0,  39.0};

	static constexpr std::array<double, 11> T_TABLE_     = { 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000};
	static constexpr std::array<double, 11> sigma_TABLE_ = {0.031, 0.27, 2.05, 6.06, 12.0, 19.6,  29.6,  41.1,  54.1,  67.7,  81.5};

	template <typename Array, typename Function>
	static Array transform(const Array& array, const Function& function) {
		Array result{};
		std::transform(array.begin(), array.end(), result.begin(), function);
		return result;
	}

	template <typename Array>
	static Array logarithm(const Array& array) {
		return transform(array, [](double x) { return std::log(x); });
	}

	const std::array<double, 11> log_T_table_ = logarithm(T_TABLE_);
	const std::array<double, 11> log_sigma_table_ = logarithm(sigma_TABLE_);

	/*
	 * Interpolation
	 */
	template <typename Array>
	static double linear_interpolation(double x, const Array& xs, const Array& ys) {
		const auto lower = std::lower_bound(xs.begin(), xs.end(), x);
		const auto i = static_cast<size_t>(std::distance(xs.begin(), lower)) - static_cast<size_t>(lower != xs.begin());
		return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]);
	}

	static double T0(double i) {
		return linear_interpolation(i, I_TABLE_, T0_TABLE_);
	}

	static double m(double i) {
		return linear_interpolation(i, I_TABLE_, m_TABLE_);
	}

	double sigma(double t) const {
		return std::exp(linear_interpolation(std::log(t), log_T_table_, log_sigma_table_));
	}

	static double T(double z, double T0, double Tw, double m) {
		return T0 + (Tw - T0) * std::pow(z, m);
	}

	template <typename Function>
	static double trapezoidal(const Function& f, double x0, double xn, size_t n) {
		assert(x0 < xn && n > 1);

		const double h = (xn - x0) / static_cast<double>(n);
		double s = (f(x0) + f(xn)) / 2.0;
		for (size_t i = 1; i < n; ++i) {
			s += f(x0 + static_cast<double>(i)*h);
		}
		s *= h;

		return s;
	}

	double Rp(double i) const {
		const auto integral = trapezoidal([&](double z) { return z * sigma(T(z, T0(i), Tw_, m(i))); }, 0, 1, 512);
		return Le_ / (2 * M_PI * R_ * R_ * integral);
	}

	static double Ucp(double I, double _Rp) {
		return I * _Rp;
	}

	double f(double I, double Uc) const {
		I = std::abs(I);
		return (Uc - (Rk_ + Rp(I)) * I) / Lk_;
	}

	double g(double I) const {
		return -std::abs(I) / Ck_;
	}

	std::pair<double, double> find_next_IUc(double dt, double In, double Ucn) const {
		In = std::abs(In);

		const double k1 = f(In, Ucn);
		const double m1 = g(In);

		const double k2 = f(In + dt * k1 / 2.0, Ucn + dt * m1 / 2.0);
		const double m2 = g(In + dt * k1 / 2.0);

		const double k3 = f(In + dt * k2 / 2.0, Ucn + dt * m2 / 2.0);
		const double m3 = g(In + dt * k2 / 2.0);

		const double k4 = f(In + dt * k3 / 2.0, Ucn + dt * m3 / 2.0);
		const double m4 = g(In + dt * k3 / 2.0);

		const double Inp1 = In + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
		const double Ucnp1 = Ucn + dt * (m1 + 2*m2 + 2*m3 + m4) / 6.0;

		return {std::abs(Inp1), Ucnp1};
	}
};

Dependency solve(const Parameters& parameters) {
	return Solver(parameters).solve();
}
