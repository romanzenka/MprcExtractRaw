#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>

// Setup for extracting multiple chromatograms from the MS1 spectra
class ChromatogramExtractor
{
	// Input parameters
private:
	// File to save this output to
	std::string outputFileName;

	// a colon delimited list of precursor m/z values
	std::string precursorMZsForXIC;

	// number of decimals for rounding the m/z values
	int massRounding;

	// mass tolerance for looking up the precursor m/z
	// values in the spectral data
	double ppmMassTol;

private:
	// Internal state being built spectrum by spectrum
	struct ExtractorState;
	ExtractorState *state;

public:
	// Typedefs for the vector of doubles
	typedef std::vector<double> TDoubleVec;
	typedef TDoubleVec::const_iterator TDoubleVecIter;

	ChromatogramExtractor(std::string outputFileName, std::string &precursorMZsForXIC, int massRounding, double ppmMassTol);
	~ChromatogramExtractor();

	// Record a new spectrum
	void addSpectrum(int scanNumber, const TDoubleVec *mzs, const TDoubleVec *intensities);

	// Save the extracted chromatogram data to a file
	void save();

private:
	// Initialize the state. Called by the constructor.
	void initialize();
};

