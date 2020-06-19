#include "model.hpp"
#include <cmath>

Model::Model(const Parameters &parameters) :
		Fmax_(parameters.Fmax),
		Tmax_(parameters.Tmax),
		nu_(1.0 / parameters.nu)
{
	if (parameters.table_h) {
		table_h();
	}

	if (parameters.table_tau) {
		table_tau();
	}

	if (parameters.plot_T) {
		plot_T();
	}

	if (parameters.plot_c) {
		plot_c();
	}

	if (parameters.plot_impulse) {
		plot_impulse();
	}
}

void Model::table_h() {
	is_only_first_ = false;
	eps_run_ = 1e-1;
	eps_it_ = 1e-3;

	tau_ = 0.1;
	double step_h[] = {1, 0.1, 0.01, 0.001};

	for (int i = 0; i < 4; ++i) {
		temp = QVector<QVector<double>>();
		h_ = step_h[i];

		iterations();

		QVector<double> row;
		const int index = i ? static_cast<int>(1 / h_ - 1) : 1;
		for (auto t: temp) {
			row.append(t[index]);
		}

		testH.append(row);
	}
}

void Model::table_tau() {
	eps_run_ = 1e-2;
	eps_it_ = 1e-6;

	h_ = 0.01;
	double step_tau[] = {1, 0.1, 0.01, 0.001};
	is_only_first_ = true;

	for (int i = 0; i < 4; ++i) {
		temp = QVector<QVector<double>>();
		tau_ = step_tau[i];

		iterations();
		QVector<double> row;

		const int index = i ? static_cast<int>(1 / tau_ - 1) : 1;
		for (auto t : temp[index]) {
			row.append(t);
		}

		testTau.append(row);
	}
}

void Model::plot_T() {
	eps_run_ = 1e-2;
	eps_it_ = 1e-6;
	h_ = 0.01;
	tau_ = Tmax_ / 1000.0;
	is_only_first_ = false;
	iterations();
	testXn = temp;
}

void Model::plot_c() {
	tau_ = Tmax_ / 1000.0;
	h_ = 0.01;

	is_only_first_ = false;
	double a[] = {2.049, 5, 10, 15};
	double b[] = {0.000564, 0.001, 0.01, 0.1};
	eps_run_ = 1e-3;
	eps_it_ = 1e-6;

	for (int i = 0; i < 4; ++i) {
		temp = QVector<QVector<double>>();
		a2_ = a[i];
		b2_ = b[i];

		iterations();

		QVector<double> row;
		for (auto t : temp) {
			row.append(t[0]);
		}

		testAB.append(row);
	}
}

void Model::plot_impulse() {
	is_only_first_ = false;
	is_impulse_ = true;
	eps_run_ = 1e-1;
	eps_it_ = 1e-4;
	h_ = 0.01;
	tau_ = Tmax_ / 1000.0;

	temp = QVector<QVector<double>>();
	iterations();
	for (auto t : temp) {
		testImpulse.append(t[0]);
	}
}


constexpr double Model::alpha(double x) {
	constexpr double _d = (alphaN_ * l_) / (alphaN_ - alpha0_);
	constexpr double _c = -alpha0_ * _d;
	return _c / (x - _d);
}

double Model::k(double T) {
	return a1_ * (b1_ + c1_ * std::pow(T, m1_));
}

double Model::c(double T) const {
	return a2_ + b2_ * std::pow(T, m2_) - c2_ / (T * T);
}

constexpr double Model::p(double x) {
	return 2 * alpha(x) / R_;
}

constexpr double Model::f(double x) {
	return 2 * alpha(x) * T0_ / R_;
}

double Model::chi_nph(double T) const {
	return (k(T) + k(T + tau_)) / 2;
}

double Model::chi_nmh(double T) const {
	return (k(T) + k(T - tau_)) / 2;
}

double Model::A(double T) {
	return chi_nmh(T) * tau_ / h_;
}

double Model::B(double T, double x) {
	return A(T) + D(T) + c(T) * h_ + p(x) * h_ * tau_;
}

double Model::D(double T) {
	return chi_nph(T) * tau_ / h_;
}

double Model::F(double T, double x) {
	return f(x) * h_ * tau_ + c(T) * T * h_;
}

double Model::c_nph(double T) const {
	return (c(T) + c(T + tau_)) / 2;
}

double Model::c_nmh(double T) const {
	return (c(T) + c(T - tau_)) / 2;
}

double Model::p_nph(double x) const {
	return (p(x) + p(x + h_)) / 2;
}

double Model::p_nmh(double x) const {
	return (p(x) + p(x - h_)) / 2;
}

double Model::f_nph(double x) const {
	return (f(x) + f(x + h_)) / 2;
}

double Model::f_nmh(double x) const {
	return (f(x) + f(x - h_)) / 2;
}

void Model::iterations() {
	QVector<double> tZero;
	int n = int(l_ / h_) + 1;
	for (int i = 0; i < n; ++i) {
		tZero.append(T0_);
	}

	temp.append(tZero);

	double t = tau_;

	do {
		QVector<double> prev;
		QVector<double> curr = temp.last();

		do {
			prev = curr;
			curr = sweep(prev, t);
		} while (!endRunTrought(prev, curr));

		t += tau_;
		temp.append(curr);

		if (is_only_first_ && t >= 1)
			break;
	} while (!endIterations(t));
}

std::tuple<double, double, double> Model::left_border(const QVector<double> &prev, double t) const {
	const double p_0 = p(0);
	const double p_h = p_nph(0);
	const double f_0 = f(0);
	const double f_h = f_nph(0);

	const double c_0   = c(prev[0]);
	const double c_h   = c_nph(prev[0]);
	const double chi_h = chi_nph(prev[0]);

	const double K0 = h_ * c_h / 8 + h_ * c_0 / 4 + tau_ * chi_h / h_ + tau_ * h_ * p_h / 8 + tau_ * h_ * p_0 / 4;
	const double M0 = h_ * c_h / 8 - tau_ * chi_h / h_ + tau_ * h_ * p_h / 8;
	const double P0 = h_ * c_h / 8 * (prev[0] + prev[1]) + h_ * c_0 / 4 * prev[0] + F(t) * tau_ + tau_ * h_ / 4 * (f_h + f_0);

	return {K0, M0, P0};
}

std::tuple<double, double, double> Model::right_border(const QVector<double> &prev) const {
	const double p_N   = p(l_);
	const double p_Nmh = p_nmh(l_);
	const double f_N   = f(l_);
	const double f_Nmh = (f_N + f(l_ - h_)) / 2;

	const int N = prev.size() - 1;

	const double c_N     = c(prev[N]);
	const double c_Nmh   = c_nmh(prev[N]);
	const double chi_Nmh = chi_nmh(prev[N]);

	const double KN = h_ * c_N / 4 + h_ * c_Nmh / 8 + chi_Nmh * tau_ / h_ + alphaN_ * tau_ + p_N * tau_ * h_ / 4 + p_Nmh * tau_ * h_ / 8;
	const double MN = h_ * c_Nmh / 8 - chi_Nmh * tau_ / h_ + p_Nmh * tau_ * h_ / 8;
	const double PN = h_ * c_N * prev[N] / 4 + h_ * c_Nmh * (prev[N] + prev[N - 1]) / 8 + alphaN_ * T0_ * tau_ + (f_N + f_Nmh) * tau_ * h_ / 4;

	return {KN, MN, PN};
}


QVector<double> Model::sweep(const QVector<double> &prev, const double t) {
	const auto [K0, M0, P0] = left_border(prev, t);
	const auto [KN, MN, PN] = right_border(prev);

	// forward sweep
	QVector<double> xi = {0, -M0 / K0};
	QVector<double> eta = {0, P0 / K0};
	for (auto [x, n] = std::tuple{h_, 1}; x + h_ < l_; x += h_, ++n) {
		const double xiN = xi.last();
		const double etaN = eta.last();
		const double det = (B(prev[n], x) - A(prev[n]) * xiN);
		xi.append(D(prev[n]) / det);
		eta.append((F(prev[n], x) + A(prev[n]) * etaN) / det);
	}

	// backward substitution
	QVector<double> y(xi.count());
	y[y.count() - 1] = (PN - MN * eta.last()) / (KN + MN * xi.last());
	for (int i = y.count() - 2; i >= 0; --i) {
		y[i] = xi[i + 1] * y[i + 1] + eta[i + 1];
	}

	return y;
}

bool Model::endIterations(const double t) {
	if (is_impulse_) {
		return t > 300;
	}

	int last = temp.count() - 1;
	for (int i = 0; i < temp[last].count(); ++i) {
		if (std::fabs((temp[last][i] - temp[last - 1][i]) / temp[last][i]) > eps_it_)
			return false;
	}

	return true;
}

bool Model::endRunTrought(
	const QVector<double> &prev,
	const QVector<double> &current
) {
	double max = std::fabs((current[0] - prev[0]) / current[0]);

	for (int i = 1; i < std::min(current.count(), prev.count()); ++i) {
		double e = std::fabs((current[i] - prev[i]) / current[i]);

		if (e > max)
			max = e;
	}

	return max < eps_run_;
}

double Model::F(double t) const {
	if (is_impulse_) {
		t = std::round(std::fmod(t, nu_) * 100) / 100.0;
		if (nu_ <= 0.01) t = 10;
	}

	return (Fmax_ / Tmax_) * t * std::exp(1 - t / Tmax_);
}
