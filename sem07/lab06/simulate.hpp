#pragma once

struct SimulateParams {
	struct Client {
		double t, dt;
		int N;
	} client;

	struct Terminal {
		double t, dt;
		int L;
	} terminal;

	struct Window {
		double t, dt;
		double p_return;
		int L;
	} window1, window2, window3;
};

struct SimulateResults {
	struct Device {
		int n_dropped;
		double p_dropped;
	} terminal, window1, window2, window3;
};

SimulateResults simulateEvent(const SimulateParams& params);
