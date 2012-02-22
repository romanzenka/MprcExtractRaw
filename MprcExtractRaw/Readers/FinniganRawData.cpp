// Written by Navdeep Jaitly for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// E-mail: navdeep.jaitly@pnl.gov
// Website: http://ncrr.pnl.gov/software
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#include ".\finniganrawdata.h"
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits>
#include <ATLComTime.h>
#include "..\Utilities\Interpolation.h"

#ifdef XCALIBUR_INSTALLED
namespace Engine 
{
	namespace Readers
	{
		FinniganRawData::~FinniganRawData(void)
		{			
			if (marr_data_block != 0)
			{
				delete [] marr_data_block ;
				delete [] marr_temp_data_block ;
			} 
			if (marr_rawfileName != NULL)
			{
				delete [] marr_rawfileName ; 
				marr_rawfileName = 0 ; 
			}
			// Let the COM object go
			if(m_xraw2_class != NULL)
				m_xraw2_class->Release();
		};

		FinniganRawData::FinniganRawData(void)
		{
			m_xraw2_class = NULL ;
			marr_data_block = NULL ; 
			marr_temp_data_block = NULL ; 
			mint_last_scan_size = 0 ; 
			marr_rawfileName = new char[512] ; 
			// Obtain an instance of the Finnigan file reader ActiveX control
			long nRet = GetXRawFileInstance();
			if(nRet)
				std::cerr << "Unable to get instance of XRawFile.ocx" << std::endl ;
		};

		void FinniganRawData::GetScanDescription(int scan, char *description)
		{
			_bstr_t bstr_filter ;
			m_xraw2_class->GetFilterForScanNum((long)scan, &bstr_filter.GetBSTR()) ;  
			strcpy(description,(char*)bstr_filter);		
		}


		const char* FinniganRawData::GetFileName()
		{
			return marr_rawfileName ; 
		}

		int FinniganRawData::GetScanSize()
		{
			return mint_last_scan_size ; 
		}

		int FinniganRawData::GetNumScans()	
		{ 
			return (int)mlong_num_spectra ; 
		}
		int FinniganRawData::FirstSpectraNumber() 
		{ 
			return (int)mlong_spectra_num_first ; 
		}
		int FinniganRawData::LastSpectraNumber() 
		{ 
			return (int)mlong_spectra_num_last ; 
		}
		int FinniganRawData::Open(char *raw_file_name)
		{
			strcpy(marr_rawfileName, raw_file_name) ; 
			
			_bstr_t bstr = marr_rawfileName ;
			BSTR bstr_type = bstr.GetBSTR() ;

			HRESULT res = m_xraw2_class->Open(bstr_type) ; 
			if(res != S_OK)
			{
				char message[512] ; 
				strcpy(message, "Unable to open XCalibur file: ") ; 
				strcat(message, marr_rawfileName) ; 
				throw message ;
			}

			// Get the number of spectra
			res = m_xraw2_class->SetCurrentController(0,1) ; 
			long nRet = m_xraw2_class->GetNumSpectra (&mlong_num_spectra );
			nRet = m_xraw2_class->GetFirstSpectrumNumber (&mlong_spectra_num_first) ;
			nRet = m_xraw2_class->GetLastSpectrumNumber (&mlong_spectra_num_last) ;
			
			
			
			if( nRet )
			{
				std::cerr << "Unable to get number of spectra from " << marr_rawfileName <<std::endl;
				return 1;
			}


			return 0;
		}


		void FinniganRawData::Close()
		{
			m_xraw2_class->Close();
		}

		void FinniganRawData::Load(char *file_n)
		{
			Open(file_n);
		}

		double FinniganRawData::GetSignalRange(int scan_num)
		{
			if (scan_num == mint_last_scan_num)
				return mdbl_signal_range ; 

			int lastWholeNumber = 0;
			double signal_range = 0 ; 
			VARIANT varMassList;
			VariantInit(&varMassList);
			VARIANT varPeakFlags;
			VariantInit(&varPeakFlags);
			long nArraySize = 0;

			_bstr_t bstr = "" ;
			BSTR bstr_type = bstr.GetBSTR() ; 
			double peak_width ;

			long scanN = scan_num;

			HRESULT res = m_xraw2_class->SetCurrentController(0,1) ;
			long nRet = m_xraw2_class->GetMassListFromScanNum (&scanN, 
						bstr_type,	// no filter
						0,			// no cutoff
						0,			// no cutoff
						0,			// all peaks returned
						FALSE,			// do not centroid
						&peak_width,
						&varMassList,		// mass list data
						&varPeakFlags,		// peak flags data
						&nArraySize );		// size of mass list array

			if( nArraySize )
			{
				// Get a pointer to the SafeArray
				SAFEARRAY FAR* psa = varMassList.parray;

				DataPeak* pDataPeaks = NULL;
				SafeArrayAccessData( psa, (void**)(&pDataPeaks) );

				double min_intensity = DBL_MAX ; 
				double max_intensity = DBL_MIN ; 

				for( long j=0; j<nArraySize; j++ )
				{
					double intensity = pDataPeaks[j].dIntensity ; 
					if (intensity > max_intensity) 
						max_intensity = intensity ; 
					if (intensity < min_intensity) 
						min_intensity = intensity ; 
				}
			
				signal_range = (max_intensity - min_intensity) ; 
			}

			if( varMassList.vt != VT_EMPTY )
			{
				SAFEARRAY FAR* psa = varMassList.parray;
				varMassList.parray = NULL;

				// Delete the SafeArray
				SafeArrayDestroy( psa );
			}

			if(varPeakFlags.vt != VT_EMPTY )
			{
				SAFEARRAY FAR* psa = varPeakFlags.parray;
				varPeakFlags.parray = NULL;

				// Delete the SafeArray
				SafeArrayDestroy( psa );
			}

			return signal_range ; 

		}

		double FinniganRawData::GetScanTime(int scan_num)
		{
			if (scan_num == mint_last_scan_num)
				return mdbl_last_scan_time ; 
			double start_time ; 

			m_xraw2_class->RTFromScanNum((long)scan_num, &start_time) ; 
			return start_time ; 
		}

	
		bool FinniganRawData::GetRawData(std::vector<double> *mzs, std::vector<double> *intensities, int scan_num, int num_points)
		{
			// Finnigan data is already truncated. Dont mess with it. 
			return GetRawData(mzs, intensities, scan_num) ; 
		}
	
		int FinniganRawData::GetParentScan(int scan_num)
		{			
			int msN_level = GetMSLevel(scan_num) ;
			int level = msN_level;

			int i = 0;

			while (level >= msN_level)
			{
				level = GetMSLevel(scan_num - i);
				i++;
			}
			i--;
			return(scan_num - i);
		}

		bool FinniganRawData::IsProfileScan(int scan_num)
		{

			long ms_level = 0;
			m_xraw2_class->IsProfileScanForScanNum(scan_num, &ms_level) ;
			if (ms_level == 1)
				return true ;
			else
				return false ;
		}

		bool FinniganRawData::IsFTScan(int scan_num)
		{	
			char ch_filter [512] ; 					
			_bstr_t bstr_filter ;
			m_xraw2_class->GetFilterForScanNum((long)scan_num, &bstr_filter.GetBSTR()) ;  
			strcpy(ch_filter,(char*)bstr_filter);		
			if (_strnicmp(ch_filter, "ft",2)== 0)
			{
				return true ; 
			}
			return false ;			
		}

		double FinniganRawData::GetAGCAccumulationTime(int scan_num)
		{			
			variant_t var_value ; 
			m_xraw2_class->GetTrailerExtraValueForScanNum((long)scan_num, "Ion Injection Time (ms):", &var_value.GetVARIANT());
			double time = 0.0 ; 
			time = (double) var_value.fltVal ; 
			return time ; 
		}

		double FinniganRawData::GetTICForScan(int scan_num)
		{
			long num_packets ; 
			double start_time ; 
			double low_mass ;
			double high_mass ; 
			double tic ; 
			double base_peak ; 
			double base_intensity ; 
			double frequency ; 
			long unif_time ; 
			long num_channels ; 

			m_xraw2_class->GetScanHeaderInfoForScanNum(scan_num, &num_packets, &start_time, &low_mass, &high_mass, &tic, &base_peak,
					&base_intensity, &num_channels, &unif_time, &frequency) ;  
			
			return tic ; 
		}

		// Fills the buf with dissociation type. The buffer must be at least 10 characters long (+trailing zero)
		void FinniganRawData::GetDissociationType(int scan_num, char *buf) {
			char ch_filter [512] ; 					
			_bstr_t bstr_filter ;
			m_xraw2_class->GetFilterForScanNum((long)scan_num, &bstr_filter.GetBSTR()) ;  			
			strcpy(ch_filter,(char*)bstr_filter);		
			int bufPos=0;
			char* in = strstr(ch_filter, "@");
			
			if(in==NULL) {
				*buf=0;
				return;
			}

			in++;
			char* out = buf;
			while(*in && !isdigit(*in) && out-buf<10) {
				*out++=*in++;
			}
			*out=0;
		}

		int FinniganRawData::GetMSLevel(int scan_num)
		{
			int ms_level = 1;

			//gets the filter string
			char ch_filter [512] ; 					
			_bstr_t bstr_filter ;
			m_xraw2_class->GetFilterForScanNum((long)scan_num, &bstr_filter.GetBSTR()) ;  
			strcpy(ch_filter,(char*)bstr_filter);		
			
			//search for 'ms'
			for ( int chNum = 0; chNum < 512; chNum++)
			{				
				if (ch_filter[chNum] == 'm')
				{
					if(ch_filter[chNum+1] == 's')
					{
						chNum = chNum+2;
						char ch = ch_filter[chNum] ;
						char ch1 = ch_filter[chNum+1] ; 
						int ms = (int) ch ; 
						switch(ch)
						{
							case '2': ms_level = 2;
									  break;
							case '3': ms_level = 3;
									  break;
							case '4': ms_level = 4 ; 
									  break ;
							case '5': ms_level = 5 ; 
									  break ; 
							case '6': ms_level = 6 ; 
									  break ; 
							case '7': ms_level = 7 ; 
									  break ; 
							case '8': ms_level = 8 ; 
									  break ; 
							case '9': ms_level = 9 ; 
									  break ; 
							case '1': ms_level = 0 ; 
									break ;
						/*	case '10': ms_level  = 10 ; 
									   break ; 
							case '11': ms_level = 11 ; 
									   break ; 
							case '12': ms_level = 12 ; 
									   break ;
							case '13': ms_level  = 13 ; 
									   break ; 
							case '14': ms_level  = 14 ; 
									   break ; 
							case '15': ms_level  = 15 ; 
									   break ; */
							case ' ': ms_level = 1;
									  break;
							default : ms_level = 1;
									  break;
						}

						if (ms_level == 0)
						{
							switch(ch1)
							{
								case '0': ms_level  = 10 ; 
									   break ; 
								case '1': ms_level = 11 ; 
										break ; 
								case '2': ms_level = 12 ; 
										break ;
								case '3': ms_level  = 13 ; 
										break ; 
								case '4': ms_level  = 14 ; 
										break ; 
								case '5': ms_level  = 15 ; 
										break ; 
								case '6': ms_level  = 16 ; 
										break ;
								case '7': ms_level  = 17 ; 
										break ;
								case ' ':ms_level = 1 ; 
										break ; 
								default :ms_level = 1;
										break;
								
							}
						}
						return ms_level;
					}											
					
				}
			}
			return ms_level ; 
			
		
		}
		
		bool FinniganRawData::IsMSScan (int scan_num)
		{
			//Returns true if the scan is a MS-level scan

			int ms_level = GetMSLevel(scan_num);
			if (ms_level == 1)
				return true;
			else
				return false;
		}	
		
		double FinniganRawData::GetMonoMZFromHeader(int scan_num)
		{
			double mono_mz = 0 ; 
			variant_t var_value ; 
			m_xraw2_class->GetTrailerExtraValueForScanNum((long)scan_num, "Monoisotopic M/Z:", &var_value.GetVARIANT());
			mono_mz = var_value.dblVal ; 
			return mono_mz ; 				
		}

		short FinniganRawData::GetMonoChargeFromHeader(int scan_num) 
		{
			short cs = 0 ; 
			variant_t var_value ; 
			m_xraw2_class->GetTrailerExtraValueForScanNum((long)scan_num, "Charge State:", &var_value.GetVARIANT());
			cs = (short)var_value.iVal ; 
			return cs ;
		}


		
		double FinniganRawData::GetParentMz(int scan_num)
		{
			//Return the parent m/z of the particular msN scan
			double parent_mz = 0;
			
			//gets the filter string
			char ch_filter [512] ; 		
			char ch_mz[16];
			_bstr_t bstr_filter ;
			m_xraw2_class->GetFilterForScanNum((long)scan_num, &bstr_filter.GetBSTR()) ;  
			strcpy(ch_filter,(char*)bstr_filter);		

 			int ms_level = GetMSLevel(scan_num);
			
			int parent_count = 0;
			
			if (ms_level == 2)
			{
				for ( int chNum = 0; chNum < 512; chNum++)
				{					
					if (ch_filter[chNum] == '2')
					{
						chNum++;
						int mzIndex = 0;
						while(ch_filter[chNum] != '@')
						{
							ch_mz[mzIndex] = ch_filter[chNum];
							chNum++;
							mzIndex++;
						}
						break;
					}
				}
				
			}
			else
			{
				int chNum;
				for (chNum = 0; chNum < 512; chNum++)
				{
					if (ch_filter[chNum] == '@')
					{
						parent_count++;
						if (parent_count <= (ms_level - 1))
							break;
					}
				}
				
				while(ch_filter[chNum] != ' ')
						chNum ++;
				
				int mzIndex = 0;
				while(ch_filter[chNum] != '@')
				{
					ch_mz[mzIndex] = ch_filter[chNum];
					chNum++;
					mzIndex++;
				}
				
			}

			parent_mz = atof(ch_mz);		
			
			return parent_mz;

		}
		
		bool FinniganRawData::GetRawData(std::vector<double> *mzs, std::vector<double> *intensities, int scan_num)
		{
			mint_last_scan_num = scan_num ; 
			int lastWholeNumber = 0;

			VARIANT varMassList;
			VariantInit(&varMassList);
			VARIANT varPeakFlags;
			VariantInit(&varPeakFlags);
			long nArraySize = 0;

			_bstr_t bstr = "" ;
			BSTR bstr_type = bstr.GetBSTR() ; 
			double peak_width ;

			long scanN = scan_num;

			HRESULT res = m_xraw2_class->SetCurrentController(0,1) ;
			long nRet = m_xraw2_class->GetMassListFromScanNum (&scanN, 
						bstr_type,	// no filter
						0,			// no cutoff
						0,			// no cutoff
						0,			// all peaks returned
						FALSE,			// do not centroid
						&peak_width,
						&varMassList,		// mass list data
						&varPeakFlags,		// peak flags data
						&nArraySize );		// size of mass list array

			long num_packets ; 
			double start_time ; 
			double low_mass ;
			double high_mass ; 
			double tic ; 
			double base_peak ; 
			double base_intensity ; 
			double frequency ; 
			long unif_time ; 
			long num_channels ; 

			m_xraw2_class->GetScanHeaderInfoForScanNum(scanN, &num_packets, &start_time, &low_mass, &high_mass, &tic, &base_peak,
				&base_intensity, &num_channels, &unif_time, &frequency) ;  
			mdbl_last_scan_time = start_time ; 
			if( nArraySize )
			{
				// Get a pointer to the SafeArray
				SAFEARRAY FAR* psa = varMassList.parray;

				DataPeak* pDataPeaks = NULL;
				SafeArrayAccessData( psa, (void**)(&pDataPeaks) );

				intensities->clear();
				mzs->clear();
				if ( nArraySize < (int)intensities->capacity())
				{
					intensities->reserve(nArraySize) ;
					mzs->reserve(nArraySize);
				}

				double min_intensity = DBL_MAX ; 
				double max_intensity = DBL_MIN ; 

				for( long j=0; j<nArraySize; j++ )
				{
					if (pDataPeaks[j].dMass > high_mass)
					{
						break ; 
					}
					double intensity = pDataPeaks[j].dIntensity ; 
					if (intensity > max_intensity) 
						max_intensity = intensity ; 
					if (intensity < min_intensity) 
						min_intensity = intensity ; 

					mzs->push_back(pDataPeaks[j].dMass) ;
					intensities->push_back(intensity) ;
				}
			
				mdbl_signal_range = (max_intensity - min_intensity) ; 

				// Release the data handle
				SafeArrayUnaccessData( psa );
			}

			if( varMassList.vt != VT_EMPTY )
			{
				SAFEARRAY FAR* psa = varMassList.parray;
				varMassList.parray = NULL;

				// Delete the SafeArray
				SafeArrayDestroy( psa );
			}

			if(varPeakFlags.vt != VT_EMPTY )
			{
				SAFEARRAY FAR* psa = varPeakFlags.parray;
				varPeakFlags.parray = NULL;

				// Delete the SafeArray
				SafeArrayDestroy( psa );
			}

			mint_last_scan_size = (int) mzs->size() ; 
			if (nArraySize == 0)
				return false ; 
			return true ; 
		}

		

		void FinniganRawData::GetTicFromFile(std::vector<double> *intensities, std::vector<double> *scan_times, bool base_peak_tic)
		{

			long num_packets ; 
			double start_time ; 
			double low_mass ;
			double high_mass ; 
			double tic ; 
			double base_peak ; 
			double base_intensity ; 
			double frequency ; 
			long unif_time ; 
			long num_channels ; 

			for (long scan_num = mlong_spectra_num_first ; scan_num <= mlong_spectra_num_last ; scan_num++)
			{
				m_xraw2_class->GetScanHeaderInfoForScanNum(scan_num, &num_packets, &start_time, &low_mass, &high_mass, &tic, &base_peak,
					&base_intensity, &num_channels, &unif_time, &frequency) ;  
				if (base_peak_tic)
					intensities->push_back(base_intensity) ; 
				else
					intensities->push_back(tic) ; 
				scan_times->push_back(start_time) ; 
			}

		}

		// Obtain information about the scan
		void FinniganRawData::GetScanInfo(long scan_num, double *tic, double *retentionTime, 
				double *lowMass, double *highMass,
				int *childScans, 
				double *ionInjectionTimeMs, double *cycleTimeSeconds, double *elapsedTimeSeconds, 
				double *timeToNextScanSeconds, 
				bool *lockMassFound, double *lockMassShift, double *conI, double *conA, double *conB, double *conC, double *conD, double *conE,
				double *sourceCurrent,
				double *vacuumIonGauge,
				double *vacuumConvectronGauge,
				double *ftVacuumPenningGauge,
				double *ftVacuumPiraniGauge1,
				double *ionMultiplier1,
				double *ionMultiplier2,
				double *ftCeMeasureVoltage,
				double *ftAnalyzerTemp)
		{
			using namespace std;

			long numPackets ; // How many mass/intensity pairs
			double basePeakMass ; // Mass of the base peak
			double basePeakIntensity ; // Intensity of the base peak
			long numChannels ; // Number of channels acquired for the scan
			long uniformTime ; // Indicating whether the sampling time increment for the current controller is uniform (?)
			double frequency ; // Sampling frequency for the current controller			

			m_xraw2_class->GetScanHeaderInfoForScanNum(scan_num, &numPackets, retentionTime, lowMass, highMass, tic, &basePeakMass,
				&basePeakIntensity, &numChannels, &uniformTime, &frequency) ;			

			// Extract delta information by parsing subsequent scans
			*timeToNextScanSeconds = 0;
			*cycleTimeSeconds = 0;
			*elapsedTimeSeconds = 0;
			*childScans = 0;
			if(GetMSLevel(scan_num)==1) {
				int nextScan = GetNextScanNum(scan_num);
				while(nextScan <= mlong_spectra_num_last) {
					int msLevel = GetMSLevel(nextScan);
					double nextRT = GetScanTime(nextScan);

					if(msLevel==1) {
						*cycleTimeSeconds = (nextRT-*retentionTime)*60.0;
						break;
					} else {
						*childScans = *childScans+1;
					}

					nextScan = GetNextScanNum(nextScan);
				}
			} 
			
			if(scan_num+1<=mlong_spectra_num_last) {
				double nextRT = GetScanTime(scan_num+1);
				*timeToNextScanSeconds = (nextRT-*retentionTime)*60.0;
			}

			this->GetScanHeaderData(scan_num, ionInjectionTimeMs, elapsedTimeSeconds, lockMassFound, lockMassShift, conI, conA, conB, conC, conD, conE);
			this->GetStatusLogData(scan_num,
				sourceCurrent,
				vacuumIonGauge,
				vacuumConvectronGauge,
				ftVacuumPenningGauge,
				ftVacuumPiraniGauge1,
				ionMultiplier1,
				ionMultiplier2,
				ftCeMeasureVoltage,
				ftAnalyzerTemp);
		}

		void FinniganRawData::GetScanHeaderData(long scan_num, double *ionInjectionTimeMs, double *elapsedTimeSeconds, bool *lockMassFound, double *lockMassShift, double *conI, double *conA, double *conB, double *conC, double *conD, double *conE) {
			// Try to obtain additional values from the trailer data
			*ionInjectionTimeMs = 0;
			*lockMassFound = false;
			*lockMassShift = 0.0;		
			*conI = *conA = *conB = *conC = *conD = *conE = 0.0;			

			VARIANT varLabels;
			VariantInit(&varLabels);
			VARIANT varValues;
			VariantInit(&varValues);
			long nArraySize = 0;						
			long nRet = m_xraw2_class->GetTrailerExtraForScanNum(scan_num,
				&varLabels,
				&varValues,
				&nArraySize);

			// Get a pointer to the SafeArray
			SAFEARRAY FAR* psaLabels = varLabels.parray;
			varLabels.parray = NULL;
			SAFEARRAY FAR* psaValues = varValues.parray;
			varValues.parray = NULL;
			BSTR* pbstrLabels = NULL;
			BSTR* pbstrValues = NULL;
			if( FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels) ) ) )
			{
				SafeArrayUnaccessData( psaLabels );
				SafeArrayDestroy( psaLabels );
				throw "Failed to access labels array for scan trailer info";
			}
			if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
			{
				SafeArrayUnaccessData( psaLabels );
				SafeArrayDestroy( psaLabels );
				SafeArrayUnaccessData( psaValues );
				SafeArrayDestroy( psaValues );
				throw "Failed to access values array for scan trailer info";
			}
			for( long i=0; i<nArraySize; i++ )
			{
				std::string sLabel = _bstr_t(pbstrLabels[i]);
				std::string sData = _bstr_t(pbstrValues[i]);
								
				if(sLabel=="Ion Injection Time (ms):") {
					*ionInjectionTimeMs = atof(sData.c_str());
				} else if(sLabel=="Elapsed Scan Time (sec):") {
					*elapsedTimeSeconds = atof(sData.c_str());
				} else if(sLabel=="Conversion Parameter I:") {
					*conI = atof(sData.c_str());
				} else if(sLabel=="Conversion Parameter A:") {
					*conA = atof(sData.c_str());
				} else if(sLabel=="Conversion Parameter B:") {
					*conB = atof(sData.c_str());
				} else if(sLabel=="Conversion Parameter C:") {
					*conC = atof(sData.c_str());
				} else if(sLabel=="Conversion Parameter D:") {
					*conD = atof(sData.c_str());
				} else if(sLabel=="Conversion Parameter E:") {
					*conE = atof(sData.c_str());
				} else if(sLabel=="FT Analyzer Message:") {
					const char * lockPartStart = "Lock(";
					const char * lockPartEnd = ")";
					const char * unitSuffix = "ppm";
					// Try to parse out information about lockmass
					int lockStart = sData.rfind(lockPartStart);
					if(lockStart!=std::string::npos) {
						int lockEnd = sData.find(lockPartEnd, lockStart+strlen(lockPartStart));
						if(lockEnd!=std::string::npos) {
							std::string lockInfo = sData.substr(lockStart+strlen(lockPartStart), lockEnd-lockStart-strlen(lockPartStart));
							if(lockInfo.substr(lockInfo.length()-3)==unitSuffix) {								
								// We can extract lockmass info. Our string ends with +-<num>ppm
								if(lockInfo.find("NOT FOUND")==std::string::npos) {
									*lockMassFound=true;
								}
								int lmNumStart=lockInfo.rfind(' ');
								if(lmNumStart!=std::string::npos) {
									std::string ppmVal = lockInfo.substr(lmNumStart+1, lockInfo.length()-strlen(unitSuffix)-(lmNumStart+1));									
									*lockMassShift=atof(ppmVal.c_str());
								}								
							}
						}
					}
				}
			}
			// Delete the SafeArray
			SafeArrayUnaccessData( psaLabels );
			SafeArrayDestroy( psaLabels );
			SafeArrayUnaccessData( psaValues );
			SafeArrayDestroy( psaValues );			
		}

		void FinniganRawData::GetStatusLogData(long scan_num, 
			double *sourceCurrent,
			double *vacuumIonGauge,
			double *vacuumConvectronGauge,
			double *ftVacuumPenningGauge,
			double *ftVacuumPiraniGauge1,
			double *ionMultiplier1,
			double *ionMultiplier2,
			double *ftCeMeasureVoltage,
			double *ftAnalyzerTemp) {
			// Try to obtain additional values from the status log

			*sourceCurrent = 0.0;
			*vacuumIonGauge = 0.0;
			*vacuumConvectronGauge = 0.0;
			*ftVacuumPenningGauge = 0.0;
			*ftVacuumPiraniGauge1 = 0.0;
			*ionMultiplier1 = 0.0;
			*ionMultiplier2 = 0.0;
			*ftCeMeasureVoltage = 0.0;
			*ftAnalyzerTemp = 0.0;

			// -- Boilerplate
			VARIANT varLabels;
			VariantInit(&varLabels);
			VARIANT varValues;
			VariantInit(&varValues);
			// --

			long nArraySize = 0;						
			double statusLogRT = 0.0;
			long nRet = m_xraw2_class->GetStatusLogForScanNum(scan_num,
				&statusLogRT,
				&varLabels,
				&varValues,
				&nArraySize);

			// -- Boilerplate: Get a pointer to the SafeArray
			SAFEARRAY FAR* psaLabels = varLabels.parray;
			varLabels.parray = NULL;
			SAFEARRAY FAR* psaValues = varValues.parray;
			varValues.parray = NULL;
			BSTR* pbstrLabels = NULL;
			BSTR* pbstrValues = NULL;
			if( FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels) ) ) )
			{
				SafeArrayUnaccessData( psaLabels );
				SafeArrayDestroy( psaLabels );
				throw "Failed to access labels array for scan trailer info";
			}
			if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
			{
				SafeArrayUnaccessData( psaLabels );
				SafeArrayDestroy( psaLabels );
				SafeArrayUnaccessData( psaValues );
				SafeArrayDestroy( psaValues );
				throw "Failed to access values array for scan trailer info";
			}
			// --

			int section = 0;
			for( long i=0; i<nArraySize; i++ )
			{
				std::string sLabel = _bstr_t(pbstrLabels[i]);
				std::string sData = _bstr_t(pbstrValues[i]);								

				if(sLabel=="API SOURCE") {
					section = 1;
				} else if(sLabel=="VACUUM") {
					section = 2;
				} else if(sLabel=="FT VACUUM") {
					section = 3;
				} else if(sLabel=="TURBO PUMP") {
					section = 4;
				} else if(sLabel=="FT TURBO PUMP 1") {
					section = 5;
				} else if(sLabel=="FT TURBO PUMP 2") {
					section = 6;
				} else if(sLabel=="FT TURBO PUMP 3") {
					section = 7;
				} else if(sLabel=="ION OPTICS") {
					section = 8;
				} else if(sLabel=="MAIN RF") {
					section = 9;
				} else if(sLabel=="ION DETECTION SYSTEM") {
					section = 10;
				} else if(sLabel=="FT Analyzer") {
					section = 11;
				} else if(sLabel=="POWER SUPPLIES") {
					section = 12;
				} else if(sLabel=="FT POWER SUPPLIES") {
					section = 13;
				} else if(sLabel=="INSTRUMENT STATUS") {
					section = 14;
				} else if(sLabel=="SYRINGE PUMP") {
					section = 15;
				} else if(sLabel=="DIVERT VALVE") {
					section = 16;
				}

			switch(section) {
			case 1: // API SOURCE
				if(sLabel=="Source Current (uA):") {
					*sourceCurrent = atof(sData.c_str());
				}						
				break;
			case 2: // VACUUM
				if(sLabel=="Ion Gauge (E-5 Torr):") {
					*vacuumIonGauge = atof(sData.c_str());
				} else if(sLabel=="Convectron Gauge (Torr):") {
					*vacuumConvectronGauge = atof(sData.c_str());
				}
				break;
			case 3: // FT VACUUM
				if(sLabel=="FT Penning Gauge (E-10 Torr):") {
					*ftVacuumPenningGauge = atof(sData.c_str());
				} else if(sLabel=="FT Pirani Gauge 1 (Torr):") {
					*ftVacuumPiraniGauge1 = atof(sData.c_str());
				}
				break;
			case 10:  // ION DETECTION SYSTEM
				if(sLabel=="Multiplier 1 (V):") {
					*ionMultiplier1 = atof(sData.c_str());
				} else if(sLabel=="Multiplier 2 (V):") {
					*ionMultiplier2 = atof(sData.c_str());
				}
				break;
			case 11: // FT Analyzer
				if(sLabel=="FT CE Measure Voltage (V):") {
					*ftCeMeasureVoltage = atof(sData.c_str());
				} else if(sLabel.compare(0, strlen("FT Analyzer Temp. ("), "FT Analyzer Temp. (")==0) {
					*ftAnalyzerTemp = atof(sData.c_str());
				}
				break;
			}
				
			}
			
			// -- Boilerplate Delete the SafeArray
			SafeArrayUnaccessData( psaLabels );
			SafeArrayDestroy( psaLabels );
			SafeArrayUnaccessData( psaValues );
			SafeArrayDestroy( psaValues );			
			// --
		}

		void FinniganRawData::GetRawFileInfo(			
			std::string *originalFileName, 
			std::string *instrumentName, 
			std::string *instrumentSerial, 
			std::string *creationDate, 
			double *runTimeInSeconds,
			std::string *comment,			
			std::string *sampleId,
			long *numMs1,
			long *numMs2,
			long *numMs3Plus)
		{
			using namespace std;

			BSTR strBSTR = NULL;
			m_xraw2_class->GetSeqRowRawFileName(&strBSTR);
			*originalFileName = _bstr_t(strBSTR);
			
			DATE dt = NULL;
			m_xraw2_class->GetCreationDate(&dt);
			COleDateTime date(dt);
			//Convert DATE to a formatted string
			*creationDate = (LPCTSTR)date.Format("%Y-%m-%d %H:%M:%S");

			double endTimeInMinutes = 0.0;
			m_xraw2_class->GetEndTime(&endTimeInMinutes);

			*runTimeInSeconds = endTimeInMinutes*60.0;

			strBSTR = NULL;
			m_xraw2_class->GetSeqRowSampleID(&strBSTR);
			*sampleId = _bstr_t(strBSTR);
			
			strBSTR = NULL;
			m_xraw2_class->GetInstSerialNumber(&strBSTR);
			*instrumentSerial = _bstr_t(strBSTR);

			strBSTR = NULL;
			m_xraw2_class->GetInstName(&strBSTR);
			*instrumentName = _bstr_t(strBSTR);

			strBSTR = NULL;
			m_xraw2_class->GetSeqRowComment(&strBSTR);
			*comment = _bstr_t(strBSTR);

			CountSpectraByMsLevel(numMs1, numMs2, numMs3Plus);
		}

		void FinniganRawData::GetTuneMethod(std::string *tuneMethod) {
		}

		void GetInstrumentMethod(std::string *instrumentMethod) {
		}

		void GetSampleInformation(std::string *sampleInformation) {
		}

		void GetErrorLog(std::string *errorLog) {
		}

		void FinniganRawData::CountSpectraByMsLevel(long *ms1, long *ms2, long *ms3plus) {
			int first=GetFirstScanNum();
			int last=GetLastScanNum();
			for(int i=first; i<=last; i++) {
				switch(GetMSLevel(i)) {
				case 1:
					(*ms1)++;
					break;
				case 2:
					(*ms2)++;
					break;
				default:
					(*ms3plus)++;
					break;
				}
			}
		}

		int FinniganRawData::GetXRawFileInstance(void)
		{
			CoInitialize( NULL );
			CLSID clsid ; 

			HRESULT res =  CLSIDFromString(L"{5FE970B2-29C3-11D3-811D-00104B304896}", &clsid ); 	

			if (res ==REGDB_E_WRITEREGDB)
			{
				throw "Unable to instantiate Finnigan objects: XRawFile from XRawfile.dll (version 2.0.0). Please check that the following dlls from XCalibur are avaialble on your system: CFRDBResources.dll, CFRUtil.dll, ExploreDataObjects.dll, ExploreDataObjectsManaged.dll, ExploreDataObjectsps.dll, FControl2.dll, Fglobal.dll, Fileio.dll, finDB.dll, finSSClientLib.dll, Fregistry.dll, XRawfile2.dll" ;
			}
			IID riid ;

			res = IIDFromString(L"{5FE970B1-29C3-11D3-811D-00104B304896}", &riid) ; 

			if (res == E_INVALIDARG)
			{
				throw "Unable to instantiate Finnigan objects: XRawFile from XRawfile.dll (version 2.0.0). Please check that the following dlls from XCalibur are avaialble on your system: CFRDBResources.dll, CFRUtil.dll, ExploreDataObjects.dll, ExploreDataObjectsManaged.dll, ExploreDataObjectsps.dll, FControl2.dll, Fglobal.dll, Fileio.dll, finDB.dll, finSSClientLib.dll, Fregistry.dll, XRawfile2.dll" ;
			}

			res = CoCreateInstance(clsid, NULL,CLSCTX_INPROC_SERVER, riid, (void **) &m_xraw2_class);
			if(res != S_OK)
			{
				res = MyCoCreateInstance(LPCTSTR("XRawfile2.dll"), clsid, NULL, riid, (void **) &m_xraw2_class);
				
				if(res != S_OK)
				{
					throw "Unable to instantiate Finnigan objects: XRawFile from XRawfile.dll (version 2.0.0). Please check that the following dlls from XCalibur are avaialble on your system: CFRDBResources.dll, CFRUtil.dll, ExploreDataObjects.dll, ExploreDataObjectsManaged.dll, ExploreDataObjectsps.dll, FControl2.dll, Fglobal.dll, Fileio.dll, finDB.dll, finSSClientLib.dll, Fregistry.dll, XRawfile2.dll" ;
				}
			}
			
			return 0;
		}

		HRESULT __stdcall FinniganRawData::MyCoCreateInstance(LPCTSTR szDllName, IN REFCLSID rclsid, IUnknown* pUnkOuter, IN REFIID riid, OUT LPVOID FAR* ppv)
		{
		  HRESULT hr = REGDB_E_KEYMISSING;

		  HMODULE hDll = ::LoadLibrary(szDllName);
		  if (hDll == 0)
			return hr;

		  typedef HRESULT (__stdcall *pDllGetClassObject)(IN REFCLSID rclsid, IN REFIID riid, OUT LPVOID FAR* ppv);

		  pDllGetClassObject GetClassObject = (pDllGetClassObject)::GetProcAddress(hDll, "DllGetClassObject");

		  if (GetClassObject == 0)
		  {
			::FreeLibrary(hDll);
			return hr;
		  }

		  IClassFactory *pIFactory;

		  hr = GetClassObject(rclsid, IID_IClassFactory, (LPVOID *)&pIFactory);

		  if (!SUCCEEDED(hr))
			return hr;

		  hr = pIFactory->CreateInstance(pUnkOuter, riid, ppv);
		  pIFactory->Release();

		  return hr;
		}
	}
}
#endif
