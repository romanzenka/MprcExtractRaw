#include "PolymerDetection.h"
#include <algorithm> // HASA: copy

const double PolymerDetection::FRAGMENT_TOLERANCE = 0.25;
const double PolymerDetection::PROTON_MASS = 1.00727646677;
const double PolymerDetection::BIN_SIZE = 100.0;
const double PolymerDetection::MIN_INTENSITY = 5.0;

void PolymerDetection::FindPolymers(double MgfMz, int MgfZ, const std::vector<double> *inMzs, const std::vector<double> *inIntensities, double &maxOffset, double &maxSegment, double &maxScore, double &maxpValue)
{
		// Go through all possible theoretical polymer models, score them, output the scores as a table
		// A polymer model is characterized by segment_size and offset. The polymer peaks are then at
		//    offset, offset+segment_size, offset+2*segment_size, ...
		// The maximum peak mass has to be less than the reported M/Z. 
		// We test only singly and doubly-charged peaks (unless the spectrum is singly-charged)
		// For each model, we calculate cross-correlation score which is
		// sum of peak scores (peak has score 6 when it is intensity #1 within its 100 Da bin)

		static std::vector<double> mzs;
		static std::vector<int> scores;
		maxScore = 0.0;
		maxSegment = -1;
		maxOffset = -1;
		maxpValue = 1.0;
		if(inMzs->size()==0) return;

		FilterImportantPeaks(inMzs, inIntensities, BIN_SIZE, &mzs, &scores);

		double mass = MzToMass(MgfMz, MgfZ);


		double mzSpan = inMzs->at(inMzs->size()-1) - inMzs->at(0);

		for (int segment_size = MIN_SEGMENT_SIZE; segment_size <= MAX_SEGMENT_SIZE; segment_size++)
		{
			// How many possible hits can we get
			int numSegments = (int)((mzSpan / segment_size) + 1);
			std::vector<double> &pValues = *GetPValues(numSegments);
			int numpValues = pValues.size();

			for (int offset = 0; offset < segment_size; offset++)
			{
				size_t score = 0;

				for (int charge = 2; charge <= 2; charge++)
				{
					double offset_Z = MassToMz(offset, charge);
					double max_Mz = MassToMz(mass, charge);
					double segment_size_Z = segment_size / (double)charge;

					std::vector<int>::const_iterator scoresIter = scores.cbegin();
					for (std::vector<double>::const_iterator mzsIter = mzs.cbegin(); mzsIter!=mzs.cend(); mzsIter++, scoresIter++)
					{
						double mzsPeak = *mzsIter;
						if (mzsPeak < max_Mz)
						{
							int segment_num = ((int)((mzsPeak - offset_Z) / segment_size_Z + 0.5));
							double ref_peak = offset_Z + segment_size_Z * segment_num;
							if (fabs(ref_peak - mzsPeak) < FRAGMENT_TOLERANCE)
							{
								score += *scoresIter;
							}
						}						
					}
				}

				if (score >= numpValues) { score = numpValues - 1; }
				if (score > maxScore || (score == maxScore && pValues[score] < maxpValue))
				{
					maxScore = score;
					maxpValue = pValues[score];
					maxSegment = segment_size;
					maxOffset = offset;
				}
			}
		}
	}

// Retains only peaksPerInterval peaks per each intervalSize. The peak intensities
// are replaced with peak scores (peak with intensity #1 has highest score)
void PolymerDetection::FilterImportantPeaks(const std::vector<double> *mzs, const std::vector<double> *intensities, 
	double intervalSize,
	std::vector<double> *mzsOut, std::vector<int> *scoresOut)
{
	mzsOut->clear();
	scoresOut->clear();

	if(mzs->size()==0) return;

	double minMz = mzs->at(0);
	double maxMz = mzs->at(mzs->size()-1);
	int numIntervals = (int)((maxMz - minMz) / intervalSize + 1);
	int totalPeaks = NUM_PEAKS * numIntervals;

	int prevInterval = -1;
	int topIndices[NUM_PEAKS];
	double smallestIntensity = 0.0;

	int i = 0;

	while (true)
	{
		int interval;
		double peakMz;
		double peakIntensity;

		if (i==mzs->size()) 
		{
			interval = -1;
		} 
		else 
		{
			peakMz = mzs->at(i);
			peakIntensity = intensities->at(i);
			interval = (int)((peakMz - minMz) / intervalSize);
		}

		// Last entry triggers dumping the interval and exit		
		if (prevInterval != interval)
		{
			// Dump current list of top peaks
			if (prevInterval != -1)
			{
				for (int a = 0; a < NUM_PEAKS; a++)
				{
					if (topIndices[a] != -1)
					{
						mzsOut->push_back(mzs->at(topIndices[a]));
						scoresOut->push_back(NUM_PEAKS - a);
					}
				}
			}

			// Beyond the last peak in input data?
			if (i == mzs->size())
			{
				break;
			}

			// Clear the list of top peaks
			for (int a = 0; a < NUM_PEAKS; a++)
			{
				topIndices[a] = -1;
			}
			smallestIntensity = 0.0;
			prevInterval = interval;
		}

		if (peakIntensity > smallestIntensity && peakIntensity > MIN_INTENSITY)
		{
			for (int b = 0; b < NUM_PEAKS; b++)
			{
				if (topIndices[b] == -1)
				{
					topIndices[b] = i;
					smallestIntensity = peakIntensity;
					break;
				}
				if (peakIntensity > intensities->at(topIndices[b]))
				{
					memcpy(topIndices+b+1, topIndices+b, sizeof(int)*(NUM_PEAKS - b - 1));
					topIndices[b] = i;
					if (topIndices[NUM_PEAKS - 1] != -1)
					{
						smallestIntensity = intensities->at(topIndices[NUM_PEAKS - 1]);
					}
					break;
				}
			}
		}
		i++;
	}
}

PolymerDetection::TpValueMap PolymerDetection::pValueMap;

std::vector<double> *PolymerDetection::GetPValues(int numSegments) {
	TpValueMap::iterator pos = pValueMap.find(numSegments);
	std::vector<double> *result;
	if(pos==pValueMap.end()) {
		// Not found
		result = new std::vector<double>();
		CalculatePValues(*result, numSegments, NUM_PEAKS, BIN_SIZE);
		pValueMap[numSegments] = result;
	} else {
		result = pos->second;
	}	
	return result;
}

// Generates pValues for each score that can be attained for numSegments, with numPeaks retained per each binSize.
// It determines the probability that score will be at least X by random chance.
// The score of most abundant peak is numPeaks, scores go down as the peak goes down in intensity.
// Each peak is considered to occupy 1Da of the binSize.
void PolymerDetection::CalculatePValues(std::vector<double> &pValues, int numSegments, int numPeaks, double binSize)
{
	int maxScore = numSegments * numPeaks;

	pValues.resize(maxScore + 1); // pValues[i] = Probability that score > i
	std::vector<double> prevPValues;
	prevPValues.resize(maxScore+1);
	double pHit = 1.0 / binSize; // Probability that we hit a peak of a specific score
	double pMiss = (binSize - numPeaks) / binSize; // Probability we do not hit any peak

	// First segment probabilities
	prevPValues[0] = pMiss; // We get 0 if we miss
	for (int score = 1; score <= numPeaks; score++)
	{
		prevPValues[score] = pHit; // We get the score only if we hit peak of the proper rank
	}

	// Add probabilities taking into account more segments
	for (int segment = 1; segment < numSegments; segment++)
	{
		for (int score = 0; score <= maxScore; score++)
		{
			// We have no hit in this segment
			pValues[score] = pMiss * prevPValues[score];

			// We have a hit on peak that provides addScore
			for (int addScore = 1; addScore <= __min(numPeaks, score); addScore++)
			{
				pValues[score] += pHit * prevPValues[score - addScore];
			}
		}
		prevPValues.assign(pValues.begin(), pValues.end());
	}

	// Resulting pValue is a sum of all probabilities greater or equal to given score
	// (What is the probability we get at least score X by chance?)
	double sum = 0.0;
	for (int score = maxScore; score >= 0; score--)
	{
		sum += prevPValues[score];
		pValues[score] = sum;
	}
}
