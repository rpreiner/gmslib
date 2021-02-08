//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// Copyright (c) Reinhold Preiner 2014-2020 
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------

#pragma once

#include "ext/MemoryUsage.h"


class Memory
{
private:
	double mBaseline;
	double mPeakUsage;
	bool mActive;

public:
	static Memory* instance()
	{
		static Memory* instance = nullptr;
		if (!instance) 
			instance = new Memory();
		return instance;
	}
	
	void setActive(bool active)
	{
		mActive = active;
	}

	bool isActive()
	{
		return mActive;
	}

	Memory() : mActive(false)
	{
		reset(); 
	}

	void reset()
	{
		mBaseline = mPeakUsage = 0;
	}

	void setBaseline()
	{
		if (mActive)
			mBaseline = double(MemoryInfo::Usage()) / (1 << 20);
	}

	double baseline()
	{
		return mBaseline;
	}

	double record(int linenr)
	{
		if (!mActive)
			return 0;
		
		double mem = double(MemoryInfo::Usage()) / (1 << 20);
		if (mem > mPeakUsage)
		{
			mPeakUsage = mem;
			cout << "new mempeak @" << linenr << ": " << mem << " MB" << endl;
		}
		return mem;
	}

	double peakUsage()
	{
		return mPeakUsage;
	}

	double peakRelativeUsage()
	{
		return mPeakUsage - mBaseline;
	}
};


