#include "simulate.hpp"

#include <algorithm>
#include <functional>
#include <set>
#include "random.hpp"

class Processor {
public:
	virtual ~Processor() {}

	virtual void receive_request() = 0;
};

class RequestGenerator {
public:
	using Generator = std::function<double()>;

	RequestGenerator(Generator generator)
		: generator_(generator)
		, receivers_()
	{ }

	void add_receiver(Processor* receiver) {
		receivers_.insert(receiver);
	}

	double next_time_period() const {
		return generator_();
	}

	void emit_request() const {
		for (auto&& receiver : receivers_) {
			receiver->receive_request();
		}
	}

private:
	Generator generator_;
	std::set<Processor*> receivers_;
};

class RequestProcessor : public RequestGenerator, public Processor {
public:
	using Generator = RequestGenerator::Generator;

	RequestProcessor(Generator generator, double p_return)
		: RequestGenerator(generator)
		, Processor()
		, cur_queue_size_(0)
		, max_queue_size_(0)
		, n_processed_requests_(0)
		, p_return_(p_return)
		, n_return_(0)
	{ }

	~RequestProcessor() override {}

	int n_processed_requests() const { return n_processed_requests_; }
	int cur_queue_size() const { return cur_queue_size_; }
	int max_queue_size() const { return max_queue_size_; }
	int n_return() const { return n_return_; }

	void process() {
		if (cur_queue_size_ > 0) {
			++n_processed_requests_;
			--cur_queue_size_;
			emit_request();
			if (uniform_real(0.0, 1.0) < p_return_) {
				++n_return_;
				receive_request();
			}
		}
	}

	void receive_request() override {
		++cur_queue_size_;
		max_queue_size_ = std::max(cur_queue_size_, max_queue_size_);
	}

private:
	int cur_queue_size_;
	int max_queue_size_;
	int n_processed_requests_;
	const double p_return_;
	int n_return_;
};

SimulateResults simulateDt(const SimulateParams& params) {
	RequestGenerator generator([=]() { return uniform_real(params.a, params.b); });
	RequestProcessor processor([=]() { return normal(params.mu, params.sigma); }, params.p_return);
	generator.add_receiver(&processor);

	auto gen_period = generator.next_time_period();
	auto proc_period = gen_period + processor.next_time_period();
	for (double t = 0.0; processor.n_processed_requests() < params.n; t += params.dt) {
		if (gen_period <= t) {
			generator.emit_request();
			gen_period += generator.next_time_period();
		}
		if (t >= proc_period) {
			processor.process();
			if (processor.cur_queue_size() > 0) {
				proc_period += processor.next_time_period();
			} else {
				proc_period = gen_period + processor.next_time_period();
			}
		}
	}
	return {processor.n_return(), processor.max_queue_size()};
}

SimulateResults simulateEvent(const SimulateParams& params) {
	RequestGenerator generator([=]() { return uniform_real(params.a, params.b); });
	RequestProcessor processor([=]() { return normal(params.mu, params.sigma); }, params.p_return);
	generator.add_receiver(&processor);

	auto gen_period = generator.next_time_period();
	auto proc_period = gen_period + processor.next_time_period();
	while (processor.n_processed_requests() < params.n) {
		if (gen_period <= proc_period) {
			generator.emit_request();
			gen_period += generator.next_time_period();
		}
		if (gen_period >= proc_period) {
			processor.process();
			if (processor.cur_queue_size() > 0) {
				proc_period += processor.next_time_period();
			} else {
				proc_period = gen_period + processor.next_time_period();
			}
		}
	}
	return {processor.n_return(), processor.max_queue_size()};
}
