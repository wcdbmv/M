#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <QVector>

class Parameters {
public:
	double F0;
};

class Dependency {
public:
	QVector<double> x;
	QVector<QVector<double>> Tx;
	QVector<double> t;
	QVector<QVector<double>> Tt;
};

//Dependency solve(const Parameters& parameters);
std::pair<QVector<double>, QVector<std::pair<double, QVector<double>>>> solve(const Parameters& parameters);

#endif  // SOLVE_HPP_
