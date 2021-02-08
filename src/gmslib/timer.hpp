//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// Copyright (c) Reinhold Preiner 2014-2020 
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------


#pragma once

#include <chrono>

class Timer
{
	using clock = std::chrono::high_resolution_clock;

private:
	clock::time_point mStartTime;
	bool mActive = true;

public:
	Timer()
	{
		start();
	}

	void setActive(bool active)
	{
		mActive = active;
	}

	// reset timer and start counting
	void start()
	{
		if (!mActive) return;
		mStartTime = clock::now();
	}

	// returns current time relative to last start() call (in ms)
	double stop()
	{
		if (!mActive) return 0;

		return std::chrono::duration<double, std::milli>(clock::now() - mStartTime).count();
	}
};

