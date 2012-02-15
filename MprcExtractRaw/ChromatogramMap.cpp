#include "ChromatogramMap.h"

#include <stdio.h>
#include <iostream>
#include <sstream>
#include "gif/gif_lib.h"

ChromatogramMap::ChromatogramMap(int massResolution, double minMass, double maxMass)
{
	this->massResolution = massResolution;
	this->minMass = minMass;
	this->maxMass = maxMass;
	maxLogIntensity = 0.0f;
}


ChromatogramMap::~ChromatogramMap(void)
{
	for(TData::iterator it = data.begin(); it != data.end(); it++) 
	{
		delete [] (*it);
	}
	data.clear();
}

void ChromatogramMap::addSpectrum(const TMzs *inMzs, const TIntensities *inIntensities) 
{
	float *values = new float[massResolution];
	memset(values, 0, massResolution * sizeof(float));

	TMzs::const_iterator mzs=inMzs->begin();
	TIntensities::const_iterator intensities=inIntensities->begin();
	const double massRange = maxMass-minMass;
	const bool normalize = false;
	float maxSpectrumIntensity = 0.0;

	const TMzs::const_iterator mzsEnd = inMzs->end();
	while(mzs!=mzsEnd) 
	{
		double mz = *mzs;

		// Where does the value go into our array
		int bucket = (int)(((mz-minMass)/massRange)*massResolution);
		if (bucket>=0 && bucket<massResolution) 
		{
			// Intensity is log-scaled, has to be >= 1
			double intensity = *intensities;
			if(intensity<MIN_INTENSITY) 
			{ 
				intensity=MIN_INTENSITY; 
			} 
			float logIntensity = (float)log(intensity);
			if(logIntensity > maxLogIntensity) 
			{
				maxLogIntensity = logIntensity;
			}

			if(normalize) {
				if(logIntensity > maxSpectrumIntensity) 
				{
					maxSpectrumIntensity = logIntensity;
				}
				maxLogIntensity = 1.0;
			}

			// We store the maximum intensity per bucket
			if(values[bucket]<logIntensity) 
			{
				values[bucket]=logIntensity;
			}
		}

		++mzs;
		++intensities;
	}

	if(normalize) {
		for(int i=0; i<massResolution; i++) {
			values[i] = values[i]/maxSpectrumIntensity;
		}
	}

	data.push_back(values);
}

/******************************************************************************
* Handle last GIF error. Try to close the file and free all allocated memory. *
******************************************************************************/
static void HandleGifError(GifFileType *GifFile)
{
    int i = GifLastError();

    if (EGifCloseFile(GifFile) == GIF_ERROR) {
	GifLastError();
    }
    throw "Could not write out the chromatogram map";
} 

void ChromatogramMap::dumpEqualized(const std::string filename) 
{
	// Obtain a detailed histogram
	int histogramData[DETAILED_HISTOGRAM_BUCKETS];
	memset(histogramData, 0, sizeof(histogramData));	
	this->histogram(histogramData, DETAILED_HISTOGRAM_BUCKETS);

	float exportHistogram[EXPORT_HISTOGRAM_BUCKETS+1];

	// Equalize the histogram
	equalize(
		histogramData, DETAILED_HISTOGRAM_BUCKETS, data.size()*massResolution, maxLogIntensity,
		exportHistogram, EXPORT_HISTOGRAM_BUCKETS);

	// Start the .gif file
	GifFileType *gifFile;
	GifPixelType *scanLine;
	ColorMapObject *colorMap;

	if ((colorMap = MakeMapObject(256, NULL)) == NULL) {
	    throw "Failed to allocate chromatogram color map";
	}

	colorMap->Colors[0].Red=0;
	colorMap->Colors[0].Green=0;
	colorMap->Colors[0].Blue=0;
	
	for(int i=1; i<EXPORT_HISTOGRAM_BUCKETS; i++) {
		// Histogram is log-scale.
		// Convert to actual value
		double bucketCentroid = (exp(exportHistogram[i-1])+exp(exportHistogram[i]))/2.0;
		// Convert to log-10 scale - 1 is zero, 1E10 is 1
		double log10scale = log10(bucketCentroid)/10.0;
		if(log10scale<0.0) {
			log10scale=0.0;
		}
		if(log10scale>1.0) { 
			log10scale=1.0;
		}
		// Convert to 3-byte value
		long color = (long)(0xFFFFFF * log10scale);
		
		// R,G,B form together intensity, log scale, 255, 255, 255 being the maximum intensity (1E10)
		colorMap->Colors[i].Red=(unsigned char)((color & 0xFF0000)>>16);
		colorMap->Colors[i].Green=(unsigned char)((color & 0x00FF00)>>8);
		colorMap->Colors[i].Blue=(unsigned char)(color & 0x0000FF);
	}

	if ((scanLine = (GifPixelType *) malloc(sizeof(GifPixelType)*getWidth())) == NULL) {
		throw "Cannot allocate GIF scan line";
	}

	EGifSetGifVersion("89a");

	if ((gifFile = EGifOpenFileName(filename.c_str(), FALSE)) == NULL) {
		free(scanLine);
		throw "Cannot write out the chromatogram map - could not open the resulting file";
	}
	
	if (EGifPutScreenDesc(gifFile, getWidth(), getHeight(), 8, 0, colorMap) == GIF_ERROR) {
		free(scanLine);
		throw "Could not write out the chromatogram map";
	}

	if (EGifPutImageDesc(gifFile, 0, 0, getWidth(), getHeight(), FALSE, NULL) == GIF_ERROR) {
		free(scanLine);
		throw "Could not write out the chromatogram map";
	}

	// Now we can go value by value, get the resulting buckets and output a byte
	for(int row = massResolution-1; row >= 0; row--) {
		for(TData::iterator it=data.begin(); it!=data.end(); it++) 
		{
			float *spectrum = *it;

			int bucket = getBucket(*(spectrum+row), exportHistogram, EXPORT_HISTOGRAM_BUCKETS);
			char bucketChar = (char)bucket;
			
			EGifPutPixel(gifFile, bucketChar);
		}
	}

	// Create a GIF comment listing the m/z range
	std::stringstream comment;
	comment << "mz:" << minMass << "," << maxMass;
	if (EGifPutComment(gifFile, comment.str().c_str()) == GIF_ERROR) {
		free(scanLine);
		throw "Could not put comment in the file";
	}

	if (EGifCloseFile(gifFile) == GIF_ERROR) {
		free(scanLine);
		throw "Could not write out the chromatogram map";
	}
}

void ChromatogramMap::histogram(int* buckets, int numBuckets) 
{
	for(TData::iterator it = data.begin(); it != data.end(); it++) 
	{
		float *spectrum = *it;
		float *spectrumEnd = spectrum + massResolution;		
		while(spectrum < spectrumEnd) 
		{
			float intensity = *spectrum;
			int bucket = (int)((intensity/maxLogIntensity)*numBuckets);
			if(bucket>=numBuckets) 
			{
				bucket=numBuckets-1;
			}
			buckets[bucket]++;
			spectrum++;
		}
	}
}

void ChromatogramMap::equalize(int *inputBuckets, int numInputBuckets, int totalValues, float maxInputIntensity,
	float *outputBucketBoundaries, int numOutputBuckets)
{
	int consumedValues = 0;
	outputBucketBoundaries[0] = 0.0f;
	int inputBucket = 0;

	for(int i = 0; i<numOutputBuckets; i++) {
		int fill=0;
		int expectedFill = (int)(totalValues-consumedValues)/(numOutputBuckets-i);
		while (fill < expectedFill && inputBucket < numInputBuckets) {
			int current = inputBuckets[inputBucket];
			consumedValues += current;
			fill += current;
			inputBucket++;
		}
		float lastInputBucketBoundary = maxInputIntensity*inputBucket/numInputBuckets;
		outputBucketBoundaries[i+1] = lastInputBucketBoundary;
	}
}
