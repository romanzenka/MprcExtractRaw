#include "ChromatogramExtractor.h"
#include <iterator>
#include <iostream>
#include <fstream>
#include <iomanip>

// data storage for (scan#, BasePeakXIC, and ticXIC). Storage is 
// indexed by scan number
// BasePeakXIC is the maximum intensity in the given range
// ticXIC is the total intensity of all data points in the given range
typedef std::map<int, std::pair<double, double> > PeakData;

// structure to store a peptide XIC
struct PeptideChromatogram
{
	// precursor m/z of the peptide
	double precursorMZ;
	// mz window for looking up the precursor m/z in the 
	// spectral data. Window is defined as precursorMZ+/-ppmMassTol
	double minMz;
	double maxMz;
	// data storage of the XIC
	PeakData XICs;

	// default constructor
	PeptideChromatogram() {}
};

struct ChromatogramExtractor::ExtractorState {
	// Chromatograms
	std::vector<PeptideChromatogram> chromatograms;
};

double precisionRound(double value, int precision)
{
	const double adjustment = pow(10.0, precision);
	return floor(value*(adjustment)+0.5) / adjustment;
}

ChromatogramExtractor::ChromatogramExtractor(std::string outputFileName, std::string &precursorMZsForXIC, int massRounding, double ppmMassTol)
{
	this->outputFileName = outputFileName;
	this->precursorMZsForXIC = precursorMZsForXIC;
	this->massRounding = massRounding;
	this->ppmMassTol = ppmMassTol;
	initialize();
}

void ChromatogramExtractor::initialize() {
	ExtractorState *state = new ExtractorState();
	this->state = state;	

	// get the list of all precursor mzs. Each mz is separated by ':'
	// convert them into an array of doubles
	std::vector<std::string> precursorStrings;
	std::vector<double> precursors;
	char *nextToken = NULL;
	char *token = strtok_s(const_cast<char*>(this->precursorMZsForXIC.c_str()), ":", &nextToken);
	while (token != NULL) {
		if (strlen(token) >= 2) {
			precursors.push_back(atof(token));
		}
		token = strtok_s(NULL, ":", &nextToken);
	}

	// create a list of precursors for XIC
	for (std::vector<double>::iterator it = precursors.begin(); it != precursors.end();	 it++) {
		double precursorMZ = *it;
		PeptideChromatogram chrom;
		chrom.precursorMZ = precursorMZ;
		double tol = (precursorMZ / 1.0e6)*ppmMassTol;
		chrom.minMz = precursorMZ - tol;
		chrom.maxMz = precursorMZ + tol;
		state->chromatograms.push_back(chrom);
		std::cout << std::fixed << std::setprecision(massRounding) << "#Precursor Peak Info:" << std::setprecision(massRounding) << precursorMZ << "\t" << chrom.minMz << "-" << chrom.maxMz << std::endl;
	}

	std::cout << "Total number of precursors:" << state->chromatograms.size() << std::endl;
}

ChromatogramExtractor::~ChromatogramExtractor()
{
	delete this->state;
}

void ChromatogramExtractor::save() {
	// output a text file containig the summed mass chromatogram
	std::cout << "Writing peak list to " << outputFileName << std::endl;
	std::ofstream outputFileStream(outputFileName, std::ofstream::out);
	outputFileStream << "Precursor m/z\tm/z Window\tScan ID\tBasePeakXIC\tTICXIC" << std::endl;
	outputFileStream << std::setprecision(this->massRounding) << std::fixed;
	size_t totalExtractedPeaks = 0;
	for (std::vector<PeptideChromatogram>::iterator chromatogram = state->chromatograms.begin(); chromatogram != state->chromatograms.end(); chromatogram++) {
		totalExtractedPeaks += chromatogram->XICs.size();
		if (chromatogram->XICs.size() > 0) {
			for (PeakData::iterator pd = chromatogram->XICs.begin(); pd != chromatogram->XICs.end(); pd++) {
				outputFileStream << chromatogram->precursorMZ << "\t" << chromatogram->minMz << "-" << chromatogram->maxMz << "\t" << pd->first << "\t" << pd->second.first << "\t" << pd->second.second << std::endl;
			}
		}
		else {
			outputFileStream << chromatogram->precursorMZ << "\t" << chromatogram->minMz << "-" << chromatogram->maxMz << "\t" << 1 << "\t" << 0.0 << "\t" << 0.0 << std::endl;
		}
	}
	outputFileStream.flush();
	outputFileStream.close();

}

void ChromatogramExtractor::addSpectrum(int scanNumber, const TDoubleVec *mzs, const TDoubleVec *intensities) {
	// iterater over all candidate peptides
	for (std::vector<PeptideChromatogram>::iterator chromatogram = state->chromatograms.begin(); chromatogram != state->chromatograms.end(); chromatogram++) {
		TDoubleVecIter mz = mzs->begin();
		TDoubleVecIter intensity = intensities->begin();

		// for each peak, bin it according to user input and keep a total of each bin
		while (mz != mzs->end()) {			
			if (*intensity > 0) {
				double roundedMZ = *mz;
				if (chromatogram->minMz<=roundedMZ && roundedMZ<=chromatogram->maxMz) {

					chromatogram->XICs[scanNumber].first = __max(chromatogram->XICs[scanNumber].first, *intensity); // First == base peak
					chromatogram->XICs[scanNumber].second += *intensity; // Second == total ion current
				}
			}

			mz++;
			intensity++;
		}
	}
}

