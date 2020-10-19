#pragma once

#include <fstream>

class RandomNumberTable {
public:
	RandomNumberTable(std::string filename)
		: stream_(filename)
	{
		if (!stream_) {
			throw;
		}
	}

	void reopen()
	{
		stream_.clear();
		stream_.seekg(0, stream_.beg);
	}

	int operator()()
	{
		int value;
		while (!(stream_ >> value)) {
			reopen();
		}
		return value;
	}

private:
	std::ifstream stream_;
};
