#pragma once

#include "PolymerDetection.h"
#include <vector>
#include <iostream>

class Test
{
public:
	Test(void);
	~Test(void);

	void TestPolymer() {
		double mzsArr[] = { 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0,
							200.0, 300.0, 400.0 };
		double intensitiesArr[] = { 3.0, 4.0, 4.5, 4.7, 5.0, 105.0, 106.0, 107.0,
			200.0, 300.0, 400.0 };
		std::vector<double> mzs(mzsArr, mzsArr+(sizeof(mzsArr)/sizeof(double)));
		std::vector<double> intensities(intensitiesArr, intensitiesArr+(sizeof(intensitiesArr)/sizeof(double)));

		std::vector<double> mzsOut;
		std::vector<int> scoresOut;
		PolymerDetection::FilterImportantPeaks(&mzs, &intensities, 100.0, &mzsOut, &scoresOut);
		for(int i=0; i<mzsOut.size(); i++) {
			std::cout << mzsOut[i] << '\t' << scoresOut[i] << '\n';
		}
	}
};
