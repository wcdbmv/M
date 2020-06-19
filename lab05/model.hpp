#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <QVector>

struct Parameters {
	double Fmax;
	double Tmax;
	double nu;

	bool table_h;
	bool table_tau;
	bool plot_T;
	bool plot_c;
	bool plot_impulse;
};

class Model
{
public:
	explicit Model(const Parameters& parameters);
	void table_h();
	void table_tau();
	void plot_T();
	void plot_c();
	void plot_impulse();
	void iterations();

public:
	double h_;
	double tau_;

	QVector<QVector<double>> temp;

	QVector<QVector<double>> testH;
	QVector<QVector<double>> testTau;

	QVector<QVector<double>> testXn;
	QVector<QVector<double>> testAB;
	QVector<double> testImpulse;

private:
	static constexpr double alpha(double x);
	static double k(double T);
	double c(double T) const;
	static constexpr double p(double x);
	static constexpr double f(double x);
	double chi_nph(double T) const;
	double chi_nmh(double T) const;
	double A(double T);
	double B(double T, double x);
	double D(double T);
	double F(double T, double x);
	double c_nph(double T) const;
	double c_nmh(double T) const;
	double p_nph(double x) const;
	double p_nmh(double x) const;
	double f_nph(double x) const;
	double f_nmh(double x) const;
	bool endIterations(double t);
	std::tuple<double, double, double> left_border(const QVector<double> &prev, double t) const;
	std::tuple<double, double, double> right_border(const QVector<double> &prev) const;
	QVector<double> sweep(const QVector<double> &prev, const double t);
	bool endRunTrought(const QVector<double> &prev, const QVector<double> &current);
	double F(double t) const;

private:
	static constexpr double eps_ = 1e-4;
	double eps_run_;
	double eps_it_;

	static constexpr double a1_ = 0.0134;
	static constexpr double b1_ = 1;
	static constexpr double c1_ = 0.000435;
	static constexpr double m1_ = 1;

	double a2_ = 2.049;
	double b2_ = 0.000563;
	static constexpr double c2_ = 52800;
	static constexpr double m2_ = 1;

	static constexpr double alpha0_ = 0.05;
	static constexpr double alphaN_ = 0.01;
	static constexpr double l_= 10;
	static constexpr double T0_= 300;
	static constexpr double R_= 0.5;

	const double Fmax_;
	const double Tmax_;

	bool is_only_first_ = false;
	bool is_impulse_ = false;

	const double nu_;
};

#endif  // MODEL_HPP_
