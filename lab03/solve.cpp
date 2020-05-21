#include "solve.hpp"

class Solver {
public:
	explicit Solver(const Parameters& parameters) :
			k0_    (parameters.k0),
			kN_    (parameters.kN),
			alpha0_(parameters.alpha0),
			alphaN_(parameters.alphaN),
			F0_    (parameters.F0),
			b_     (l_ * kN_ / (kN_ - k0_)),
			a_     (-k0_ / b_),
			d_     (l_ * alphaN_ / (alphaN_ - alpha0_)),
			c_     (-alpha0_ * d_) {}

	Dependency solve() const {
		Container a(1), b(1), c(1), d(1), xs;
		for (double x = 0.0; x <= l_; x += h_) {
			a.push_back(x_nmh(x) / h_);
			c.push_back(x_nph(x) / h_);
			b.push_back(a.back() + c.back() + p(x) * h_);
			d.push_back(f(x) * h_);
			xs.push_back(x);
		}

		// compiler will optimize calculations
		const double k0 = x_nph(0) + h_ * h_ * (p(0) + p(h_)) / 16 + h_ * h_ * p(0) / 4;
		const double m0 = h_ * h_ * (p(0) + p(h_)) / 16 - x_nph(0);
		const double p0 = h_ * F0_ + h_ * h_ * ((f(0) + f(h_)) / 2 + f(0)) / 4;

		const double kN = -x_nmh(l_) / h_ - alphaN_ - h_ * (p(l_) + p(l_ - h_)) / 16 - h_ * p(l_) / 4;
		const double mN = x_nmh(l_) / h_ - h_ * (p(l_) + p(l_ - h_)) / 16;
		const double pN = -(alphaN_ * T0_ + h_ * ((f(l_) + f(l_ - h_)) / 2 + f(l_)) / 4);

		// forward sweep
		Container xi(a.size() + 1);
		Container eta(a.size() + 1);
		xi[1]  = -m0 / k0;
		eta[1] =  p0 / k0;
		for (int i = 1; i < a.size(); ++i) {
			const double det = b[i] - a[i] * xi[i];
			xi[i + 1]  = c[i] / det;
			eta[i + 1] = (a[i] * eta[i] + d[i]) / det;
		}

		// backward substitution
		Container ys(a.size());
		ys.back() = (pN - mN * eta.back()) / (kN + mN * xi.back());
		for (int i = ys.size() - 2; i >= 0; --i) {
			ys[i] = xi[i + 1] * ys[i + 1] + eta[i + 1];
		}

		return {xs, ys};
	}

private:
	const double k0_;
	const double kN_;
	const double alpha0_;
	const double alphaN_;
	const double F0_;

	static constexpr double l_  = 10;
	static constexpr double T0_ = 300;
	static constexpr double R_  = 0.5;

	const double b_;
	const double a_;
	const double d_;
	const double c_;

	static constexpr double h_ = 0.1;

	double k(double x) const {
		return a_ / (x - b_);
	}

	double alpha(double x) const {
		return c_ / (x - d_);
	}

	double p(double x) const {
		return 2 * alpha(x) / R_;
	}

	double f(double x) const {
		return 2 * alpha(x) * T0_ / R_;
	}

	double x_nph(double xn) const {
		const double k_curr = k(xn);
		const double k_next = k(xn + h_);
		return 2 * k_curr * k_next / (k_curr + k_next);
	}

	double x_nmh(double xn) const {
		const double k_curr = k(xn);
		const double k_prev = k(xn - h_);
		return 2 * k_curr * k_prev / (k_curr + k_prev);
	}
};

Dependency solve(const Parameters& parameters) {
	return Solver(parameters).solve();
}

