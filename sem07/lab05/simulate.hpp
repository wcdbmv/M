#pragma once

struct SimulateParams {
	int t, dt;
	int n;

	struct Operator {
		int t;
		int dt;
	} op1, op2, op3;

	struct Computer {
		int t;
	} comp1, comp2;
};

struct SimulateResults {
	int n_return;
	double p_return;
};

SimulateResults simulateEvent(const SimulateParams& params);
