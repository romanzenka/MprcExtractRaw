#pragma once

#include <vector>
#include <math.h>
#include <map>

class PolymerDetection
{
private:
	static const double FRAGMENT_TOLERANCE; // How far can the peak be from the theoretical mass (Da)
	static const double PROTON_MASS;
	
	// Peak filtering - spectrum is divided into bins and only the most intense peaks in the bin are retained
	static const int NUM_PEAKS = 6; // How many peaks per bin do we retain
	static const double BIN_SIZE; // Size of the bin used for filtering peaks
	static const double MIN_INTENSITY; // Minimal intensity the peak has to have to be considered a peak

	typedef std::map<int /* numSegments */, std::vector<double>* /* p-values */> TpValueMap;
	static TpValueMap pValueMap;

public:
	static double MassToMz(double Mass, int Z) {
		return (Mass+PROTON_MASS*Z)/Z;
	}

	static double MzToMass(double Mz, int Z) {
		return Mz * Z - PROTON_MASS * Z;
	}

	static void FindPolymers(
		// The M/Z and Z values from the spectrum
		double MgfMz, int MgfZ, 
		// m/z - intensity pairs describing the spectrum
		const std::vector<double> *inMzs, const std::vector<double> *inIntensities, 
		// What segment sizes should we find the best match for
		int minSegmentSize, int maxSegmentSize, 
		// Output values - what offest+segment do we reach the maxScore at, with what p-value
		double &maxOffset, double &maxSegment, double &maxScore, double &maxpValue);	

	// Retains only peaksPerInterval peaks per each intervalSize. The peak intensities
	// are replaced with peak scores (peak with intensity #1 has highest score)
	static void FilterImportantPeaks(const std::vector<double> *mzs, const std::vector<double> *intensities, 
		double intervalSize,
		std::vector<double> *mzsOut, std::vector<int> *scoresOut);
	
	// Returns cached pValue array that was previously calculated by CalculatePValues
	// Uses static pValueMap to store the previous results. This method makes this class not thread safe.
	static std::vector<double> *GetPValues(int numSegments);

	// Generates pValues for each score that can be attained for numSegments, with numPeaks retained.
	// The score of most abundant peak is numPeaks. Each peak is considered to occupy 1Da of the binSize.
	static void CalculatePValues(std::vector<double> &pValues, int numSegments, int numPeaks, double binSize);	
};
