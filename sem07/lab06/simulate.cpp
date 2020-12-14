#include "simulate.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <vector>

#include <QtDebug>

#include "random.hpp"

static double cur_time;

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
		, return_receiver_(nullptr)
		, n_generated_requests_(0)
		, next_event_time(0.0)
	{ }

	virtual ~RequestGenerator() { }

	int n_generated_requests() const {
		return n_generated_requests_;
	}

	void add_receiver(Processor* receiver) {
		receivers_.push_back(receiver);
	}

	void set_return_receiver(Processor* return_receiver) {
		return_receiver_ = return_receiver;
	}

	double generate_time() const {
		return generator_();
	}

	Processor* emit_request() {
		++n_generated_requests_;
		if (receivers_.size() == 3) {
			const auto i = static_cast<size_t>(uniform_int(0, 2));
			auto receiver = receivers_[i];
			if (receiver->receive_request()) {
				qDebug() << "Receive";
				return receiver;
			}
			qDebug() << "Not";
			return nullptr;
		}
		for (auto&& receiver : receivers_) {
			if (receiver->receive_request()) {
				return receiver;
			}
		}
		return nullptr;
	}

	Processor* return_request() {
		if (!return_receiver_) {
			qDebug() << "Return receiver is not defined!";
			return nullptr;
		}
		if (return_receiver_->receive_request()) {
			qDebug() << "Return receiver";
			return return_receiver_;
		}
		qDebug() << "Not Return receiver";
		return nullptr;
	}

private:
	Generator generator_;
	std::vector<Processor*> receivers_;
	Processor* return_receiver_;
	int n_generated_requests_;
public:
	double next_event_time;
};

class RequestProcessor : public RequestGenerator, public Processor {
public:
	using Generator = RequestGenerator::Generator;

	RequestProcessor(Generator generator, int max_queue_size = 0, double p_return = 0.0)
		: RequestGenerator(generator)
		, Processor()
		, max_queue_size_(max_queue_size)
		, n_queued_requests_(0)
		, n_processed_requests_(0)
		, n_dropped_requests_(0)
		, n_return_(0)
		, p_return_(p_return)
	{ }

	~RequestProcessor() override { }

	int n_queued_requests() const { return n_queued_requests_; }
	int n_processed_requests() const { return n_processed_requests_; }
	int n_dropped_requests() const { return n_dropped_requests_; }

	void process() {
		if (n_queued_requests_ > 0) {
			++n_processed_requests_;
			--n_queued_requests_;
			emit_request();

			if (uniform_real(0, 1) < p_return_) {
				++n_return_;
				return_request();
			}
		}
		next_event_time = n_queued_requests_ ? cur_time + generate_time() : 0.0;
	}

	bool receive_request() override {
		const auto cond = max_queue_size_ == 0 || n_queued_requests_ < max_queue_size_;
		if (cond) {
			++n_queued_requests_;
			next_event_time = n_queued_requests_ ? cur_time + generate_time() : 0.0;
		} else {
			++n_dropped_requests_;
		}
		return cond;
	}

private:
	int max_queue_size_;
	int n_queued_requests_;
	int n_processed_requests_;
	int n_dropped_requests_;
	int n_return_;
	double p_return_;
};

SimulateResults simulateEvent(const SimulateParams& params) {
	RequestGenerator client_generator([=]() { return uniform_real(params.client.t - params.client.dt, params.client.t + params.client.dt); });
	RequestProcessor terminal([=]() { return uniform_real(params.terminal.t - params.terminal.dt, params.terminal.t + params.terminal.dt); }, params.terminal.L);
	RequestProcessor window1([=]() { return uniform_real(params.window1.t - params.window1.dt, params.window1.t + params.window1.dt); }, params.window1.L);
	RequestProcessor window2([=]() { return uniform_real(params.window2.t - params.window2.dt, params.window2.t + params.window2.dt); }, params.window2.L);
	RequestProcessor window3([=]() { return uniform_real(params.window3.t - params.window3.dt, params.window3.t + params.window3.dt); }, params.window3.L);

	client_generator.add_receiver(&terminal);
	terminal.add_receiver(&window1);
	terminal.add_receiver(&window2);
	terminal.add_receiver(&window3);
	window1.set_return_receiver(&terminal);

	const std::array<RequestGenerator*, 5> devices{&client_generator, &terminal, &window1, &window2, &window3};

	client_generator.next_event_time = client_generator.generate_time();
	terminal.next_event_time = terminal.generate_time();
	while (client_generator.n_generated_requests() < params.client.N) {
		cur_time = client_generator.next_event_time;
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
				} else {
					client_generator.emit_request();
					client_generator.next_event_time = cur_time + client_generator.generate_time();
				}
			}
		}
	}

	const auto res = [](const RequestProcessor& processor) {
		const auto n_dropped_requests = processor.n_dropped_requests();
		const auto n_processed_requests = processor.n_processed_requests();
		const auto p_dropped_requests = static_cast<double>(n_dropped_requests) / (n_dropped_requests + n_processed_requests);
		return SimulateResults::Device{n_dropped_requests, p_dropped_requests};
	};

	return {res(terminal), res(window1), res(window2), res(window3)};
}
