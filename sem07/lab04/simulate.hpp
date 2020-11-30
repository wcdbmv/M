#pragma once

struct SimulateParams {
	double a;
	double b;

	double mu;
	double sigma;

	int n;
	double p_return;
	double dt;
};

struct SimulateResults {
	int n_return;
	int max_queue_size;
};

SimulateResults simulateDt(const SimulateParams& params);
SimulateResults simulateEvent(const SimulateParams& params);
