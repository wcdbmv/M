#include "simulate.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <set>
#include "random.hpp"

class Processor {
public:
	virtual ~Processor() {}

	virtual bool receive_request() = 0;
};

class RequestGenerator {
public:
	using Generator = std::function<double()>;

	RequestGenerator(Generator generator)
		: generator_(generator)
		, receivers_()
		, n_generated_requests_(0)
		, next_event_time(0.0)
	{ }

	virtual ~RequestGenerator() { }

	int n_generated_requests() const {
		return n_generated_requests_;
	}

	void add_receiver(Processor* receiver) {
		receivers_.insert(receiver);
	}

	double generate_time() const {
		return generator_();
	}

	Processor* emit_request() {
		++n_generated_requests_;
		for (auto&& receiver : receivers_) {
			if (receiver->receive_request()) {
				return receiver;
			}
		}
		return nullptr;
	}

private:
	Generator generator_;
	std::set<Processor*> receivers_;
	int n_generated_requests_;
public:
	double next_event_time;
};

class RequestProcessor : public RequestGenerator, public Processor {
public:
	using Generator = RequestGenerator::Generator;

	RequestProcessor(Generator generator, int max_queue_size = 0)
		: RequestGenerator(generator)
		, Processor()
		, max_queue_size_(max_queue_size)
		, n_queued_requests_(0)
		, n_processed_requests_(0)
	{ }

	~RequestProcessor() override { }

	int n_queued_requests() const { return n_queued_requests_; }
	int n_processed_requests() const { return n_processed_requests_; }

	void process() {
		if (n_queued_requests_ > 0) {
			++n_processed_requests_;
			--n_queued_requests_;
			emit_request();
		}
	}

	bool receive_request() override {
		const auto cond = max_queue_size_ == 0 || n_queued_requests_ < max_queue_size_;
		if (cond) {
			++n_queued_requests_;
		}
		return cond;
	}

private:
	int max_queue_size_;
	int n_queued_requests_;
	int n_processed_requests_;
};

SimulateResults simulateEvent(const SimulateParams& params) {
	RequestGenerator client_generator([=]() { return uniform_real(params.t - params.dt, params.t + params.dt); });
	RequestProcessor op1([=]() { return uniform_real(params.op1.t - params.op1.dt, params.op1.t + params.op1.dt); }, 1);
	RequestProcessor op2([=]() { return uniform_real(params.op2.t - params.op2.dt, params.op2.t + params.op2.dt); }, 1);
	RequestProcessor op3([=]() { return uniform_real(params.op3.t - params.op3.dt, params.op3.t + params.op3.dt); }, 1);
	RequestProcessor comp1([=]() { return params.comp1.t; });
	RequestProcessor comp2([=]() { return params.comp2.t; });

	client_generator.add_receiver(&op1);
	client_generator.add_receiver(&op2);
	client_generator.add_receiver(&op3);
	op1.add_receiver(&comp1);
	op2.add_receiver(&comp1);
	op3.add_receiver(&comp2);

	const std::array<RequestGenerator*, 6> devices{&client_generator, &op1, &op2, &op3, &comp1, &comp2};

	int n_return = 0;
	client_generator.next_event_time = client_generator.generate_time();
	op1.next_event_time = op1.generate_time();
	while (client_generator.n_generated_requests() < params.n) {
		double cur_time = client_generator.next_event_time;
		for (auto&& device : devices) {
			if (0.0 < device->next_event_time && device->next_event_time < cur_time) {
				cur_time = device->next_event_time;
			}
		}

		for (auto&& device : devices) {
			if (cur_time == device->next_event_time) {
				auto processor = dynamic_cast<RequestProcessor*>(device);
				if (processor) {
					processor->process();
					processor->next_event_time = processor->n_queued_requests() ? cur_time + processor->generate_time() : 0.0;
				} else {
					auto assigned_processor = client_generator.emit_request();
					if (assigned_processor) {
						auto pr = dynamic_cast<RequestProcessor*>(assigned_processor);
						pr->next_event_time = cur_time + pr->generate_time();
					} else {
						++n_return;
					}
					client_generator.next_event_time = cur_time + client_generator.generate_time();
				}
			}
		}
	}

	return {n_return, static_cast<double>(n_return) / params.n};
}
