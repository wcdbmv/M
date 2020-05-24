#include "solve.hpp"
#include <algorithm>
#include <cmath>

class Solver_ {
public:
	Solver_(double F0) : F0_(F0) {}

	std::pair<QVector<double>, QVector<std::pair<double, QVector<double>>>> count_T() {
		QVector<std::pair<double, QVector<double>>> results;
		QVector<double> x;

		QVector<double> prev_T;
		for (auto i = 0.; i <= l_ + h_; i += h_) {
			prev_T.push_back(300.);
			x.push_back(i);
		}
		results.push_back({0, prev_T});

		QVector<double> curr_T(prev_T.size());
		std::copy(prev_T.begin(), prev_T.end(), curr_T.begin());

		for (auto time = 0.; time <= 700; time += tau_) {
			QVector<double> previous(prev_T.size());
			std::copy(prev_T.begin(), prev_T.end(), previous.begin());

			while (true) {
				const auto [K0, M0, P0] = left_border(curr_T, prev_T);
				const auto [KN, MN, PN] = right_border(curr_T, prev_T);

				QVector<double> a(1);
				QVector<double> b(1);
				QVector<double> c(1);
				QVector<double> d(1);

				auto i = 1;
				for (auto x = 0.; x <= l_; x += h_, i++) {
					auto next_i = i + 1 == curr_T.size() ? i : i + 1;
					auto x_min_half = (k(curr_T[i]) + k(curr_T[i - 1])) / 2;
					auto x_plus_half = (k(curr_T[i]) + k(curr_T[next_i])) / 2;

					a.push_back(x_min_half * tau_ / h_);
					c.push_back(x_plus_half * tau_ / h_);
					b.push_back(-a[a.size() - 1] - c[c.size() -1]
							- this->c(curr_T[i]) * h_ - p(x) * h_ * tau_);
					d.push_back(-f(x) * h_ * tau_ - this->c(curr_T[i]) * prev_T[i] * h_);
				}

				QVector<double> etta(a.size() + 1);
				QVector<double> tetta(a.size() + 1);
				QVector<double> y(a.size());

				etta[1] = - M0 / K0;
				tetta[1] = P0 / K0;

				// direct stroke of the method
				for (auto i = 1; i < a.size(); i++) {
					auto determinator = b[i] + a[i] * etta[i];
					etta[i + 1] = -c[i] / determinator;
					tetta[i + 1] = (d[i] - a[i] * tetta[i]) / determinator;
				}

				// reverse stroke
				y[y.size() - 1] = (PN - KN * tetta[tetta.size() - 1])
						/ (MN + KN * etta[etta.size() - 1]);

				for (int i = y.size() - 2; i >= 0; i--) {
					y[i] = etta[i + 1] * y[i + 1] + tetta[i + 1];
				}

				curr_T = std::move(y);

				auto max = 0.;
				for (auto i = 0; i < curr_T.size(); i++) {
					auto diff = std::abs(curr_T[i] - previous[i]) / curr_T[i];
					max = std::max(max, diff);
				}

				if (max < eps_) {
					results.push_back({time, curr_T});
					break;
				}
				std::copy(curr_T.begin(), curr_T.end(), previous.begin());
			}

			auto max = 0.;
			for (auto i = 0; i < curr_T.size(); i++) {
				auto diff = std::abs(curr_T[i] - prev_T[i]) / curr_T[i];
				max = std::max(max, diff);
			}

			if (max < eps_) {
				break;
			}
			std::copy(curr_T.begin(), curr_T.end(), prev_T.begin());
		}

		return {x, results};
	}

private:
	const double F0_;

	static constexpr double l_   = 10;
	static constexpr double T0_  = 300;
	static constexpr double R_   = 0.5;
	static constexpr double tau_ = 2;

	static constexpr double a1_ = 0.0134;
	static constexpr double b1_ = 1;
	static constexpr double c1_ = 4.35e-4;
	static constexpr double m1_ = 1;

	static constexpr double a2_ = 2.049;
	static constexpr double b2_ = 0.563e-3;
	static constexpr double c2_ = 0.528e+5;
	static constexpr double m2_ = 1;

	static constexpr double alpha0_ = 0.05;
	static constexpr double alphaN_ = 0.01;

	static constexpr double h_ = 1e-2;
	static constexpr double eps_ = 1e-5;

	static double k(double T) {
		return a1_ * (b1_ + c1_ * std::pow(T, m1_));
	}

	static double c(double T) {
		return a2_ + b2_ * std::pow(T, m2_) - c2_ / (T * T);
	}

	static constexpr double alpha(double x) {
		constexpr double delta_alpha = alphaN_ - alpha0_;
		auto __c = -(alphaN_ * alpha0_ * l_) / delta_alpha;
		auto __d = (alphaN_ * l_) / delta_alpha;
		return __c / (x - __d);
	}

	static constexpr double p(double x) {
		return 2 * alpha(x) / R_;
	}

	static constexpr double f(double x) {
		return  2 * alpha(x) * T0_ / R_;
	}

	static double x_nph(double T_n, double T_np1) {
		return (k(T_n) + k(T_np1)) / 2;
	}

	static double x_nmh(double T_n, double T_nm1) {
		return (k(T_n) + k(T_nm1)) / 2;
	}

	std::tuple<double, double, double> left_border(const QVector<double> &T_curr, const QVector<double> &T_prev) {
		constexpr double p_0 = p(0);
		constexpr double p_h = (p_0 + p(h_)) / 2;
		constexpr double f_0 = f(0);
		constexpr double f_h = (f_0 + f(h_)) / 2;

		const double c_0   = c(T_curr[0]);
		const double c_h   = c((T_curr[0] + T_curr[1]) / 2);
		const double chi_h = (k(T_curr[0]) + k(T_curr[1])) / 2;

		const double K0 = h_ * c_h / 8 + h_ * c_0 / 4 + tau_ * chi_h / h_ + tau_ * h_ * p_h / 8 + tau_ * h_ * p_0 / 4;
		const double M0 = h_ * c_h / 8 - tau_ * chi_h / h_ + tau_ * h_ * p_h / 8;
		const double P0 = h_ * c_h / 8 * (T_prev[0] + T_prev[1]) + h_ * c_0 / 4 * T_prev[0] + F0_ * tau_ + tau_ * h_ / 4 * (f_h + f_0);

		return {K0, M0, P0};
	}

	std::tuple<double, double, double> right_border(const QVector<double> &T_curr, const QVector<double> &T_prev) {
		constexpr double p_N   = p(l_);
		constexpr double p_Nmh = (p_N + p(l_ - h_)) / 2;
		constexpr double f_N   = f(l_);
		constexpr double f_Nmh = (f_N + f(l_ - h_)) / 2;

		const int N = T_curr.size() - 1;

		const double c_N     = c(T_curr[N]);
		const double c_Nmh   = c((T_curr[N] + T_curr[N - 1]) / 2);
		const double chi_Nmh = (k(T_curr[N]) + k(T_curr[N - 1])) / 2;

		const double KN = h_ * c_N / 4 + h_ * c_Nmh / 8 + chi_Nmh * tau_ / h_ + alphaN_ * tau_ + p_N * tau_ * h_ / 4 + p_Nmh * tau_ * h_ / 8;
		const double MN = h_ * c_Nmh / 8 - chi_Nmh * tau_ / h_ + p_Nmh * tau_ * h_ / 8;
		const double PN = h_ * c_N * T_prev[N] / 4 + h_ * c_Nmh * (T_prev[N] + T_prev[N - 1]) / 8 + alphaN_ * T0_ * tau_ + (f_N + f_Nmh) * tau_ * h_ / 4;

		return {KN, MN, PN};
	}
};

class Solver {
public:
	explicit Solver(const Parameters& parameters) : F0_ (parameters.F0) {}

	Dependency solve() const {
		Dependency result;
		for (double x = 0.0; x <= l_ + h_; x += h_) {
			result.x.push_back(x);
		}
		QVector<double> T_prev(result.x.size(), T0_);
		QVector<double> T_curr(T_prev);
		result.t.push_back(0.0);
		result.Tt.push_back(T_prev);

		for (double t = tau_; t <= 700; t += tau_) {
			QVector<double> previous(T_prev);
			for (;;) {
				const auto [K0, M0, P0] = left_border(T_curr, T_prev);
				const auto [KN, MN, PN] = right_border(T_curr, T_prev);

				QVector<double> A(1), B(1), D(1), F(1);
				int i = 1;
				for (double x = 0.0; x <= l_; x += h_, ++i) {
					const int next_i = i + 1 == T_curr.size() ? i : i + 1;
					A.push_back(x_nmh(T_curr[i], T_curr[i - 1]) * tau_ / h_);
					D.push_back(x_nph(T_curr[i], T_curr[next_i]) * tau_ / h_);
					B.push_back(-A.back() - D.back() - c(T_curr[i]) * h_ - p(x) * h_ * tau_);
					F.push_back(-f(x) * h_ * tau_ - c(T_curr[i]) * T_prev[i] * h_);
				}

				// forward sweep
				QVector<double> xi(A.size() + 1);
				QVector<double> eta(A.size() + 1);
				xi[1]  = -M0 / K0;
				eta[1] =  P0 / K0;
				for (int i = 1; i < A.size(); ++i) {
					const double det = B[i] - A[i] * xi[i - 1];
					xi[i + 1]  = -D[i] / det;
					eta[i + 1] = (F[i] - A[i] * eta[i - 1]) / det;
				}

				// backward substitution
				T_curr = QVector<double>(A.size());
				T_curr.back() = (PN - KN * eta.back()) / (MN + KN * xi.back());
				for (int i = T_curr.size() - 2; i >= 0; --i) {
					T_curr[i] = xi[i + 1] * T_curr[i + 1] + eta[i + 1];
				}

				if (max_diff(T_curr, previous) < eps_) {
					result.t.push_back(t);
					result.Tt.push_back(T_curr);
					break;
				}

				qCopy(T_curr.begin(), T_curr.end(), previous.begin());
			}

			if (max_diff(T_curr, T_prev) < eps_) {
				break;
			}

			qCopy(T_curr.begin(), T_curr.end(), T_prev.begin());
		}

		const auto step = static_cast<int>(1.0 / h_);
		for (int x = 0; x <= 10; ++x) {
			result.Tx.push_back({});
			for (auto& Tt: result.Tt) {
				result.Tx.back().push_back(Tt[x * step]);
			}
		}

		return result;
	}

private:
	const double F0_;

	static constexpr double l_   = 10;
	static constexpr double T0_  = 300;
	static constexpr double R_   = 0.5;
	static constexpr double tau_ = 2;

	static constexpr double a1_ = 0.0134;
	static constexpr double b1_ = 1;
	static constexpr double c1_ = 4.35e-4;
	static constexpr double m1_ = 1;

	static constexpr double a2_ = 2.049;
	static constexpr double b2_ = 0.563e-3;
	static constexpr double c2_ = 0.528e+5;
	static constexpr double m2_ = 1;

	static constexpr double alpha0_ = 0.05;
	static constexpr double alphaN_ = 0.01;

	static constexpr double h_ = 1e-2;
	static constexpr double eps_ = 1e-5;

	static double k(double T) {
		return a1_ * (b1_ + c1_ * std::pow(T, m1_));
	}

	static double c(double T) {
		return a2_ + b2_ * std::pow(T, m2_) - c2_ / (T * T);
	}

	static constexpr double alpha(double x) {
		const double delta_alpha = alphaN_ - alpha0_;
		const double _c = -(alphaN_ * alpha0_ * l_) / delta_alpha;
		const double _d = (alphaN_ * l_) / delta_alpha;
		return _c / (x - _d);
	}

	static constexpr double p(double x) {
		return 2 * alpha(x) / R_;
	}

	static constexpr double f(double x) {
		return 2 * alpha(x) * T0_ / R_;
	}

	static double x_nph(double T_n, double T_np1) {
		return (k(T_n) + k(T_np1)) / 2;
	}

	static double x_nmh(double T_n, double T_nm1) {
		return (k(T_n) + k(T_nm1)) / 2;
	}

	// flexin'
	template <typename Container>
	static double max_diff(const Container& a, const Container& b) {
		return std::transform_reduce(
			a.begin(),
			a.end(),
			b.begin(),
			0.0,
			[](double acc, double cur) { return std::max(acc, cur); },
			[](double ai, double bi) { return std::abs((ai - bi) / ai); }
		);
	}

	std::tuple<double, double, double> left_border(const QVector<double> &T_curr, const QVector<double> &T_prev) const {
		const double c_0 = c(T_curr[0]);
		const double c_h = c((T_curr[0] + T_curr[1]) / 2);
		const double chi_h = (k(T_curr[0]) + k(T_curr[1])) / 2;
		constexpr double p_0_   = p(0);
		constexpr double p_h_   = (p_0_ + (p(h_))) / 2;
		constexpr double f_0_   = f(0);
		constexpr double f_h_   = (f_0_ + f(h_)) / 2;

		const double K0 = h_ * c_h / 8 + h_ * c_0 / 4 + tau_ * chi_h / h_ + tau_ * h_ * p_h_ / 8 + tau_ * h_ * p_0_ / 4;
		const double M0 = h_ * c_h / 8 - tau_ * chi_h / h_ + tau_ * h_ * p_h_ / 8;
		const double P0 = h_ * c_h / 8 * (T_prev[0] + T_prev[1]) + h_ * c_0 / 4 * T_prev[0] + F0_ * tau_ + tau_ * h_ / 4 * (f_h_ + f_0_);

		return {K0, M0, P0};
	}

	std::tuple<double, double, double> right_border(const QVector<double> &T_curr, const QVector<double> &T_prev) const {
		const int N = T_curr.size() - 1;

		const double c_N   = c(T_curr[N]);
		const double c_Nmh = c((T_curr[N] + T_curr[N - 1]) / 2);
		const double chi_Nmh = (k(T_curr[N]) + k(T_curr[N - 1])) / 2;
		constexpr double p_N_   = p(l_);
		constexpr double p_Nmh_ = (p_N_ + p(l_ - h_)) / 2;
		constexpr double f_N_   = f(l_);
		constexpr double f_Nmh_ = (f_N_ + f(l_ - h_)) / 2;

		const double KN = h_ * c_N / 4 + h_ * c_Nmh / 8 + chi_Nmh * tau_ / h_ + alphaN_ * tau_ + p_N_ * tau_ * h_ / 4 + p_Nmh_ * tau_ * h_ / 8;
		const double MN = h_ * c_Nmh / 8 - chi_Nmh * tau_ / h_ + p_Nmh_ * tau_ * h_ / 8;
		const double PN = h_ * c_N * T_prev[N] / 4 + h_ * c_Nmh * (T_prev[N] + T_prev[N - 1]) / 8 + alphaN_ * T0_ * tau_ + (f_N_ + f_Nmh_) * tau_ * h_ / 4;

		return {KN, MN, PN};
	}
};

//Dependency solve(const Parameters& parameters) {
//	return Solver(parameters).solve();
//}

std::pair<QVector<double>, QVector<std::pair<double, QVector<double>>>> solve(const Parameters& parameters) {
	return Solver_(parameters.F0).count_T();
}
