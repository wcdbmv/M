#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <QVector>

class Parameters {
public:
	double k0;
	double kN;
	double alpha0;
	double alphaN;
	double F0;
};

using Container = QVector<double>;

class Dependency {
public:
	Container x;
	Container T;
};

Dependency solve(const Parameters& parameters);

#endif  // SOLVE_HPP_
