#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <QVector>

struct Parameters {
	double R;
	double Le;
	double Lk;
	double Ck;
	double Rk;
	double Uc0;
	double I0;
};

using Container = QVector<double>;

struct Dependency {
	Container t;
	Container I;
	Container Uc;
	Container Rp;
	Container Ucp;
};

Dependency solve(const Parameters& parameters);

#endif  // SOLVE_HPP_
