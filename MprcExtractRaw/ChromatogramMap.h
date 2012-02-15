#pragma once

#include <vector>

class ChromatogramMap
{
private:
	static const int MIN_INTENSITY = 1;	
	static const int DETAILED_HISTOGRAM_BUCKETS=65536;
	static const int EXPORT_HISTOGRAM_BUCKETS=256;

	int massResolution; // How many pixels per the width of a spectrum
	double minMass; // Minimum mass to capture
	double maxMass; // Maximum mass to capture

	float maxLogIntensity; // Maximum log(peak intensity)
	
	typedef float *TDataRow;
	typedef std::vector<TDataRow> TData;
	TData data; // Condensed data of all spectra so far. log of each value is being stored. Values <1 are stored as 0

	typedef std::vector<double> TMzs;
	typedef std::vector<double> TIntensities;

public:
	ChromatogramMap(int massResolution, double minMass, double maxMass);
	~ChromatogramMap(void);

	// Add spectrum to the in-memory representation
	void addSpectrum(const TMzs *inMzs, const TIntensities *inIntensities);	

	// Output the equalized picture into a given file in .gif format
	void dumpEqualized(const std::string filename);

	// Return the resulting map width (== number of spectra)
	int getWidth() {
		return data.size();
	}

	// Return the resulting map height (== number of pixels per spectrum == mass resolution)
	int getHeight() {
		return massResolution;
	}


private:
	// Calculate a histogram of intensity values of given size.
	// The interval <0, maxIntensity> is split into numBuckets equal buckets.
	// We return how many values fall into a given bucket
	void histogram(int* buckets, int numBuckets);

	// Equalize a calculated histogram for given output amount of buckets
	// Inputs:
	// - inputBuckets, numInputBuckets - original histogram
	// - totalValues - total amount of entries in the original histogram (could be obtained by summing the entire histogram).
	// - maxInputIntensity - maximum intensity of the input histogram (minimum is known to == 0)
	// - outputBucketBoundaries - boundaries of the output histogram buckets (preallocated). There is one extra boundary at the end of the last bucket.
	// - numOutputBuckets - how many output buckets to produce
	static void equalize(int *inputBuckets, int numInputBuckets, int totalValues, float maxInputIntensity,
		float *outputBucketBoundaries, int numOutputBuckets);

	// Binary histogram search, returns the histogram bucket containing the given value.
	// bucket[n] = <bucketBoundaries[n], bucketBoundaries[n+1])
	// Inputs:
	// - bucketBoundaries define buckets (see above)
	// - numBuckets is one less than number of bucket boundaries
	// Returns:
	// - the bucket number encoded 
	inline int getBucket(float input, float *bucketBoundaries, int numBuckets) {
		int l = 0;
		int r = numBuckets;
		while(r-l>1) {
			int mid = (l+r)/2;
			if(bucketBoundaries[mid]>=input) {
				r=mid;
			} else {
				l=mid;
			}
		}
		return l;
	}

	inline static double clamp(double x, double a, double b) 
	{    
		return x < a ? a : (x > b ? b : x);
	}
};

