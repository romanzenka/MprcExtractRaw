#include "stdafx.h"

#include "Readers/ReaderFactory.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <io.h>
#include "SQLite/sqlite3.h"
#include "OutputBuffer.h"
#include "Readers/FinniganRawData.h"
#include "PolymerDetection.h"
#include "ChromatogramMap.h"

// Version of MprcExtractRaw
#define VERSION "0.3"

// Version of the ms_scans table
#define SPECTRA_VERSION "0.1"
#define PEAKS_VERSION "0.1"

// Command line options
const std::string mprc = "--mprc";
const std::string data = "--data";
const std::string mzRange = "--mzrange";
const std::string paramsFile = "--params";

// MPRC params
const std::string rawFile = "--raw";
const std::string outFile = "--out";
const std::string ms2Switch = "--ms2";

// Data params
const std::string infoFile = "--info";
const std::string spectraFile = "--spectra";
const std::string chromatogramFile = "--chromatogram";
const std::string tuneMethodFile = "--tune";
const std::string instrumentMethodFile = "--instrument";
const std::string sampleInformationFile = "--sample";
const std::string errorLogFile = "--errorlog";

// MZ range params (rawFile is shared with MPRC params)
const std::string minMz = "--min";
const std::string maxMz = "--max";
const std::string peaksFile = "--peaks";

// Settings
const int chromatogramMzBins = 1000; // How many bins per m/z range to report
const double secondPeakMinDistanceFromBase = 5; // How far can the second most abundant peak be from the base peak (Da)

// Column names for the spectra data file
const std::string ScanId = "Scan Id";
const std::string ParentMz = "Parent m/z";
const std::string TotalIonCurrent = "TIC";
const std::string RetentionTime = "RT";
const std::string MsLevel = "MS Level";
const std::string ParentScan = "Parent Scan";
const std::string ChildScans = "Child Scans";
const std::string IonInjectionTimeMs = "Ion Injection Time";
const std::string CycleTimeSeconds = "Cycle Time";
const std::string ElapsedTimeSeconds = "Elapsed Time";
const std::string DeadTimeSeconds = "Dead Time";
const std::string TimeToNextScanSeconds = "Time To Next Scan";
const std::string LockMassFound = "Lock Mass Found";
const std::string LockMassShift = "Lock Mass Shift";
const std::string ConI = "Conversion Parameter I";
const std::string ConA = "Conversion Parameter A";
const std::string ConB = "Conversion Parameter B";
const std::string ConC = "Conversion Parameter C";
const std::string ConD = "Conversion Parameter D";
const std::string ConE = "Conversion Parameter E";
const std::string DissociationType = "Dissociation Type";
const std::string PolymerSegment = "Polymer Segment Size";
const std::string PolymerOffset = "Polymer Offset";
const std::string PolymerScore = "Polymer Score";
const std::string PolymerPValue = "Polymer p-value";
const std::string BasePeakMz = "Base Peak m/z";
const std::string BasePeakIntensity = "Base Peak Intensity";
const std::string SecondPeakMz = "Second Peak m/z";
const std::string SecondPeakIntensity = "Second Peak Intensity";

// API SOURCE
const std::string SourceCurrent = "Source Current (uA)";
// VACUUM
const std::string VacuumIonGauge = "Vacuum Ion Gauge (E-5 Torr)";
const std::string VacuumConvectronGauge = "Vacuum Convectron Gauge (Torr)";
// FT VACUUM
const std::string FtVacuumPenningGauge = "FT Penning Gauge (E-10 Torr)";
const std::string FtVacuumPiraniGauge1 = "FT Pirani Gauge 1 (Torr)";
// ION DETECTION SYSTEM
const std::string IonMultiplier1 = "Multiplier 1 (V)";
const std::string IonMultiplier2 = "Multiplier 2 (V)";
// FT ANALYZER
const std::string FtCeMeasureVoltage = "FT CE Measure Voltage (V)";
const std::string FtAnalyzerTemp = "FT Analyzer Temp (C)";

// Parameter names for raw info data file
const std::string originalRawFileName = "Name";
const std::string creationDate = "Creation Date";
const std::string sampleId = "Sample Id";
const std::string numberOfMs1Spectra = "MS1 Spectra";
const std::string numberOfMs2Spectra = "MS2 Spectra";
const std::string numberOfMs3PlusSpectra = "MS3+ Spectra";
const std::string instrumentName = "Instrument Name";
const std::string instrumentSerial = "Instrument Serial";
const std::string runTimeInSeconds = "Run Time (seconds)";
const std::string comment = "Comment";

// Column names for the peaks data file
// ScanId from spectra
// RetentionTime from spectra
const std::string Mz = "Mz";
const std::string Intensity = "Intensity";

// Check Sqlite Error
void cse(int rc) {
	if( rc!=SQLITE_OK ){
		fprintf(stderr, "SQL error\n");
		throw 1;
	}	
}

void exec(sqlite3 *db, const char * command) {
	char *zErrMsg = 0;
	int rc = sqlite3_exec(db, command, NULL, 0, &zErrMsg);
	if( rc!=SQLITE_OK ){
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
		exit(1);
	}	
}

// Setup pragmas to improve performance
void pragmas(sqlite3 *db) {
	// Switch off synchronous mode (commands end even before data transfer is done)
	exec(db, "PRAGMA synchronous=OFF;");

	// Cache more data in memory. The value is number of 1K blocks, 2000 is default.
	exec(db, "PRAGMA cache_size=200000;");

	// Temp file is in memory
	exec(db, "PRAGMA temp_store=MEMORY;");
}

// Create version table in the database
void createVersionTable(sqlite3 *db) {
	exec(db, "CREATE TABLE version (version TEXT);");
	exec(db, "INSERT INTO version (version) values ('0.1');");
}

// Create the content and content_metadata tables 
void createContentTables(sqlite3 *db) {
	exec(db, "CREATE TABLE contents (id INTEGER PRIMARY KEY, table_name TEXT, table_version TEXT, generator TEXT, generator_version TEXT, date_created TEXT);");
	exec(db, "CREATE TABLE contents_meta (id INTEGER PRIMARY KEY, contents_id INTEGER REFERENCES contents(id), meta_key TEXT, meta_value TEXT);");
}

// Create spectra table and its metadata
void createSpectraTable(sqlite3 *db, const char *inputRawFileName) {
	exec(db, "CREATE TABLE scans (scan INTEGER PRIMARY KEY, rt REAL, ms_level INTEGER, parent_scan INTEGER, parent_mass INTEGER);");

	exec(db, "INSERT INTO contents (id, table_name, table_version, generator, generator_version, date_created) VALUES (1, 'scans', '" SPECTRA_VERSION "', 'MprcExtractRaw', '" VERSION "', date('now'));");
	char *zSQL = sqlite3_mprintf("INSERT INTO contents_meta (id, contents_id, meta_key, meta_value) values (1, 1, 'raw', '%q');", inputRawFileName);
	exec(db, zSQL);
	sqlite3_free(zSQL);
}

// Create peaks table and its metadata
void createPeaksTable(sqlite3 *db, const char *inputRawFileName, bool addMs2Peaks) {
	exec(db, "CREATE TABLE peaks (scan INTEGER PRIMARY KEY, peak_count INTEGER, peaks BLOB);");
	exec(db, "INSERT INTO contents (id, table_name, table_version, generator, generator_version, date_created) VALUES (2, 'peaks', '" PEAKS_VERSION "', 'MprcExtractRaw', '" VERSION "', date('now'));");
	char *zSQL = sqlite3_mprintf("INSERT INTO contents_meta (contents_id, meta_key, meta_value) values (2, 'has_ms2_peaks', '%q');", addMs2Peaks ? "true" : "false");
	exec(db, zSQL);
	sqlite3_free(zSQL);
}

// Add spectra from the raw data file
int addSpectra(sqlite3 *db, const char *inputRawFileName, const char * outputFileName, Engine::Readers::RawData * pRawData) {
	createSpectraTable(db, inputRawFileName);

	// Prepare statement
	const char *insertStatement = "INSERT INTO scans (scan, rt, ms_level, parent_scan, parent_mass) VALUES (?1, ?2, ?3, ?4, ?5);";

	sqlite3_stmt *pStatement=NULL;
	int rc = sqlite3_prepare_v2(db, insertStatement, strlen(insertStatement), &pStatement, NULL);
	if( rc!=SQLITE_OK ){
		std::cerr<<"SQL error when preparing insert statement.\n";
		exit(1);
	}	

	try {
		int percentDone = 0;
		int numScans = pRawData->GetNumScans();
		std::cout << "\tTotal scans: " << numScans << "\n";

		for(int i=pRawData->GetFirstScanNum(); i<=pRawData->GetLastScanNum(); i=pRawData->GetNextScanNum(i)) {
			int scanNumber = i;
			double retentionTime = pRawData->GetScanTime(scanNumber);
			int msLevel = pRawData->GetMSLevel(i); // MS level 1 indicates a survey (FTMS) scan			
			int parentScan = (msLevel>1) ? pRawData->GetParentScan(i) : 0;
			double parentMass = (msLevel>1) ? pRawData->GetParentMz(i) : 0.0;

			cse(sqlite3_bind_int(pStatement, 1, scanNumber));
			cse(sqlite3_bind_double(pStatement, 2, retentionTime));
			cse(sqlite3_bind_int(pStatement, 3, msLevel));
			cse(sqlite3_bind_int(pStatement, 4, parentScan));
			cse(sqlite3_bind_double(pStatement, 5, parentMass));

			rc = sqlite3_step(pStatement);
			if(rc!=SQLITE_DONE) {
				throw 1;
			}
			cse(sqlite3_reset(pStatement));
			int newPercent = 100 * i / numScans;
			if(newPercent!=percentDone) {
				percentDone=newPercent;
				std::cout << "Adding scan " << i << " (" << percentDone << "%)\n";
			}
		}
	} catch(...) {
		std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << outputFileName << "\n";
		return 1;
	}
	if(pStatement!=NULL) {
		cse(sqlite3_finalize(pStatement));
	}
	return 0;
}


// Add scans from the raw data file
// addMs2 enables/disables extraction of data from MS2 scans
int addMsScans(sqlite3 *db, const char *inputRawFileName, const char * outputFileName, Engine::Readers::RawData * pRawData, bool addMs2Peaks) {
	createPeaksTable(db, inputRawFileName, addMs2Peaks);

	// Prepare statement
	const char *insertStatement = "INSERT INTO peaks (scan, peak_count, peaks) VALUES (?1, ?2, ?3);";

	sqlite3_stmt *pStatement=NULL;
	int rc = sqlite3_prepare_v2(db, insertStatement, strlen(insertStatement), &pStatement, NULL);
	if( rc!=SQLITE_OK ){
		fprintf(stderr, "SQL error when preparing insert statement.\n");
		exit(1);
	}	

	OutputBuffer buf;

	try {
		int percentDone = 0;
		int numScans = pRawData->GetNumScans();

		std::vector<double> mzs;
		std::vector<double> intensities;

		for(int i=pRawData->GetFirstScanNum(); i<=pRawData->GetLastScanNum(); i=pRawData->GetNextScanNum(i)) {
			int msLevel = pRawData->GetMSLevel(i);
			// MS level 1 indicates a survey (FTMS) scan
			bool addThisScan = (msLevel==1) || (msLevel==2 && addMs2Peaks);
			if(addThisScan) {			
				pRawData->GetRawData(&mzs, &intensities, i);

				int scanNumber = i;
				int numPeaks = mzs.size();
				int numPeaksStored = 0;

				// About 20% speedup - allocating what we need in advance
				buf.ensureBufferHasSpace(numPeaks*2*sizeof(TJavaFloat));

				// OUTPUT: m/z, intensity float pairs
				std::vector<double>::iterator mzIt = mzs.begin();
				std::vector<double>::iterator intIt = intensities.begin();			
				bool wasZero = false;
				bool zeroStored = false;
				TJavaFloat previousMass = 0.0;
				while(mzIt!=mzs.end()) {
					TJavaFloat newIntensity = (TJavaFloat)*intIt;
					TJavaFloat newMass = (TJavaFloat)*mzIt;
					bool isZero = newIntensity==(TJavaFloat)0.0;

					// Skip long sequences of zeros. Retain only beginning and trailing zero for graph plotting
					if(!wasZero) {
						// Store current value
						buf.addJavaFloat(newMass);
						buf.addJavaFloat(newIntensity);
						numPeaksStored++;
						if(isZero) { // !wasZero, isZero
							wasZero=true;
							zeroStored=true;
						}
					} else { // wasZero 
						if(!isZero) {
							wasZero=false;
							if(!zeroStored) {
								// Store the previous zero
								buf.addJavaFloat(previousMass);
								buf.addJavaFloat((TJavaFloat)0.0);
								numPeaksStored++;
							}
							// Store current value
							buf.addJavaFloat(newMass);
							buf.addJavaFloat(newIntensity);
							numPeaksStored++;
						} else { // wasZero && isZero							
							// The current zero is not getting stored
							zeroStored=false;
							previousMass = newMass;	
							// Keep going, do not store anything
						}
					}
					mzIt++;
					intIt++;
				}

				// Store the last zero (if it did not get stored) to terminate our graph neatly
				if(!zeroStored) {
					buf.addJavaFloat(previousMass);
					buf.addJavaFloat((TJavaFloat)0.0);
					numPeaksStored++;
				}

				cse(sqlite3_bind_int(pStatement, 1, i));
				cse(sqlite3_bind_int(pStatement, 2, numPeaksStored));
				cse(sqlite3_bind_blob(pStatement, 3, buf.get(), buf.usedSize(), SQLITE_TRANSIENT));

				rc = sqlite3_step(pStatement);
				if(rc!=SQLITE_DONE) {
					throw 1;
				}
				cse(sqlite3_reset(pStatement));

				buf.clear();
			}		
			int newPercent = 100 * i / numScans;
			if(newPercent!=percentDone) {
				percentDone=newPercent;
				std::cout << "Adding peaks for scan " << i << " (" << percentDone << "%)\n";
			}
		}
	} catch(...) {
		std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << outputFileName << "\n";
		return 1;
	}
	if(pStatement!=NULL) {
		cse(sqlite3_finalize(pStatement));
	}
	return 0;
}

void getPolymerScore(Engine::Readers::FinniganRawData *fRawData, int scan_num, std::vector<double> *mzs, std::vector<double> *intensities, 
	double *segment, double *offset, double *score, double *pValue) {
		fRawData->GetRawData(mzs, intensities, scan_num);
		double mz = fRawData->GetMonoMZFromHeader(scan_num);
		int charge = fRawData->GetMonoChargeFromHeader(scan_num);
		PolymerDetection::FindPolymers(mz, charge, mzs, intensities, *offset, *segment, *score, *pValue);
}

void getBasePeak(int scan_num, std::vector<double> *mzs, std::vector<double> *intensities, 
	double *basePeakMz, double *basePeakIntensity, double *secondPeakMz, double *secondPeakIntensity, double minDistanceSecondFromBaseDa) { 
	*basePeakMz = 0.0;
	*basePeakIntensity = 0.0;		
	std::vector<double>::const_iterator mzsIter = mzs->cbegin();
	std::vector<double>::const_iterator intIter = intensities->cbegin();
	for (; mzsIter!=mzs->cend(); mzsIter++, intIter++) {
		if(*intIter > *basePeakIntensity) {
			*basePeakIntensity = *intIter;
			*basePeakMz = *mzsIter;
		}
	}
	*secondPeakMz = 0.0;
	*secondPeakIntensity = 0.0;
	mzsIter = mzs->cbegin();
	intIter = intensities->cbegin();
	for (;mzsIter!=mzs->cend(); mzsIter++, intIter++) {
		if(*intIter > *secondPeakIntensity && fabs(*mzsIter-*basePeakMz)>minDistanceSecondFromBaseDa) {
			*secondPeakIntensity = *intIter;
			*secondPeakMz = *mzsIter;
		}
	}
}


void extractInfoFile(Engine::Readers::FinniganRawData * fRawData, const char *infoFileName) {
	using namespace std;
	cout << "Extracting raw file information to file " << infoFileName << "." << endl;

	std::ofstream infoOutputStream;
	infoOutputStream.exceptions(std::ofstream::failbit);
	infoOutputStream.open(infoFileName);

	std::string rawFileName;
	std::string creatDate;
	std::string instName;
	std::string instSerial;
	std::string comment;
	std::string sampleId;			
	long numberOfMs1Spectra = 0;
	long numberOfMs2Spectra = 0;
	long numberOfMs3PlusSpectra = 0;
	double runTime = 0.0;

	fRawData->GetRawFileInfo(&rawFileName, &instName, &instSerial, &creatDate, &runTime, &comment, &sampleId, &numberOfMs1Spectra, &numberOfMs2Spectra, &numberOfMs3PlusSpectra);

	infoOutputStream << ::originalRawFileName << '\t' << rawFileName << endl;
	infoOutputStream << ::numberOfMs1Spectra << '\t' << numberOfMs1Spectra << endl;
	infoOutputStream << ::numberOfMs2Spectra << '\t' << numberOfMs2Spectra << endl;
	infoOutputStream << ::numberOfMs3PlusSpectra << '\t' << numberOfMs3PlusSpectra << endl;
	infoOutputStream << ::instrumentName << '\t' << instName << endl;
	infoOutputStream << ::instrumentSerial << '\t' << instSerial << endl;
	infoOutputStream << ::creationDate << '\t' << creatDate << endl;
	infoOutputStream << ::runTimeInSeconds << '\t' << runTime << endl;
	infoOutputStream << ::comment << '\t' << comment << endl;
	infoOutputStream << ::sampleId << '\t' << sampleId << endl;

	infoOutputStream.close();
}

// While we are reading spectrum by spectrum data, we can produce two files - spectrum information and chromatogram.
// Extraction of these two is thus bundled into one method for efficiency (it is very slow to list all spectra)
void extractPerSpectrumData(Engine::Readers::FinniganRawData *fRawData, std::string spectraFileName, std::string chromatogramMapFileName) {
	using namespace std;
	std::vector<double> mzs;
	std::vector<double> intensities;

	ofstream spectraOutputStream;
	if(!spectraFileName.empty()) {
		cout << "Extracting spectra information to file " << spectraFileName << "." << endl;

		spectraOutputStream.exceptions(std::ofstream::failbit);
		spectraOutputStream.open(spectraFileName);
		spectraOutputStream << std::setprecision(12);

		spectraOutputStream 
			<< ScanId << '\t' 
			<< ParentMz << '\t'
			<< TotalIonCurrent << '\t' 
			<< RetentionTime << '\t' 
			<< MsLevel << '\t' 
			<< ParentScan << '\t'
			<< ChildScans << '\t'			
			<< IonInjectionTimeMs << '\t'
			<< CycleTimeSeconds << '\t'
			<< ElapsedTimeSeconds << '\t'
			<< DeadTimeSeconds << '\t'
			<< TimeToNextScanSeconds << '\t'			
			<< LockMassFound << '\t'
			<< LockMassShift << '\t'
			<< ConI << '\t'
			<< ConA << '\t'
			<< ConB << '\t'
			<< ConC << '\t'
			<< ConD << '\t'
			<< ConE << '\t'
			<< DissociationType << '\t'

			<< PolymerSegment << '\t'
			<< PolymerOffset << '\t'
			<< PolymerScore << '\t'
			<< PolymerPValue << '\t'

			<< BasePeakMz  << '\t'
			<< BasePeakIntensity  << '\t'
			<< SecondPeakMz  << '\t'
			<< SecondPeakIntensity  << '\t'

			<< SourceCurrent << '\t'
			<< VacuumIonGauge << '\t'
			<< VacuumConvectronGauge << '\t'
			<< FtVacuumPenningGauge << '\t'
			<< FtVacuumPiraniGauge1 << '\t'
			<< IonMultiplier1 << '\t'
			<< IonMultiplier2 << '\t'
			<< FtCeMeasureVoltage << '\t'
			<< FtAnalyzerTemp
			<< endl;		
	}

	ChromatogramMap *map = NULL;

	int firstScan = fRawData->GetFirstScanNum();
	int totalScans = fRawData->GetLastScanNum()-firstScan+1;
	int scan_num = firstScan;
	int percentOutput = 0;
	bool firstSpectrum=true;
	while(scan_num <= fRawData->GetLastScanNum()) {

		double parentMz;
		double tic;
		double retentionTime;
		int msLevel;
		int parentScan;
		double lowMass;
		double highMass;
		int childScans;
		double ionInjectionTimeMs;
		double cycleTimeSeconds;
		double elapsedTimeSeconds;
		double deadTimeSeconds;
		double timeToNextScanSeconds;
		bool lockMassFound;
		double lockMassShift;
		double conI, conA, conB, conC, conD, conE;

		double polymerSegment=-1.0;
		double polymerOffset=-1.0;
		double polymerScore=0.0;
		double polymerPValue=1.0;

		double basePeakMz=0.0;
		double basePeakIntensity=0.0;
		double secondPeakMz=0.0;
		double secondPeakIntensity=0.0;

		double sourceCurrent;
		double vacuumIonGauge;
		double vacuumConvectronGauge;
		double ftVacuumPenningGauge;
		double ftVacuumPiraniGauge1;
		double ionMultiplier1;
		double ionMultiplier2;
		double ftCeMeasureVoltage;
		double ftAnalyzerTemp;

		int newPercent=100*(scan_num-firstScan)/totalScans;
		if(newPercent>percentOutput) {
			cout << newPercent << "%\n";
			percentOutput = newPercent;
		}

		fRawData->GetScanInfo(scan_num, 
			&tic, &retentionTime, &lowMass, &highMass, &childScans, &ionInjectionTimeMs, 
			&cycleTimeSeconds, &elapsedTimeSeconds,
			&timeToNextScanSeconds, 
			&lockMassFound, &lockMassShift,
			&conI, &conA, &conB, &conC, &conD, &conE,
			&sourceCurrent, &vacuumIonGauge, &vacuumConvectronGauge, &ftVacuumPenningGauge, &ftVacuumPiraniGauge1, &ionMultiplier1, &ionMultiplier2, &ftCeMeasureVoltage, &ftAnalyzerTemp);

		msLevel = fRawData->GetMSLevel(scan_num);

		parentMz = msLevel!=1 ? fRawData->GetParentMz(scan_num) : 0;

		if(msLevel==1) {
			deadTimeSeconds = cycleTimeSeconds-elapsedTimeSeconds;
			parentScan=0;
			if(map==NULL && !chromatogramMapFileName.empty()) {
				cout << "Extracting chromatogram bitmap to file " << chromatogramMapFileName << "." << endl;
				map = new ChromatogramMap(chromatogramMzBins, lowMass, highMass);
			}
			if(map!=NULL) {
				fRawData->GetRawData(&mzs, &intensities, scan_num);
				map->addSpectrum(&mzs, &intensities);
			}
		} else {
			if(!spectraFileName.empty()) {
				deadTimeSeconds = timeToNextScanSeconds - elapsedTimeSeconds;
				parentScan = fRawData->GetParentScan(scan_num);
				getPolymerScore(fRawData, scan_num, &mzs, &intensities, 
					&polymerSegment, &polymerOffset, &polymerScore, &polymerPValue);
				getBasePeak(scan_num, &mzs, &intensities, 
					&basePeakMz, &basePeakIntensity, &secondPeakMz, &secondPeakIntensity,
					secondPeakMinDistanceFromBase);
			}
		}			

		if(!spectraFileName.empty()) {
			spectraOutputStream 
				<< scan_num << '\t';

			if(parentMz!=0)
				spectraOutputStream << parentMz << '\t';
			else
				spectraOutputStream << '\t';

			spectraOutputStream 
				<< tic << '\t' 
				<< retentionTime << '\t' 				
				<< msLevel << '\t';

			if(parentScan!=0)	
				spectraOutputStream << parentScan << '\t';
			else
				spectraOutputStream << '\t';

			spectraOutputStream 
				<< childScans << '\t';

			if(firstSpectrum && (ionInjectionTimeMs == (int)ionInjectionTimeMs)) {
				// We need to make sure that we mark the first ionInjectionTime as a double, so
				// Spotfire does not get confused if many injection times look like integers.
				// Output the integer (e.g. 300) as 300.0
				spectraOutputStream << ionInjectionTimeMs << ".0\t";
			} else  {
				spectraOutputStream << ionInjectionTimeMs << '\t';
			}

			char dissociationType[11];
			fRawData->GetDissociationType(scan_num, dissociationType);

			spectraOutputStream 
				<< cycleTimeSeconds << '\t'
				<< elapsedTimeSeconds << '\t'
				<< deadTimeSeconds << '\t'				
				<< timeToNextScanSeconds << '\t'
				<< lockMassFound << '\t'
				<< lockMassShift << '\t'
				<< conI << '\t'
				<< conA << '\t'
				<< conB << '\t'
				<< conC << '\t'
				<< conD << '\t'
				<< conE << '\t'
				<< dissociationType;

			// Polymers
			spectraOutputStream << '\t'
				<< polymerSegment << '\t'
				<< polymerOffset << '\t'
				<< polymerScore << '\t';
			if(firstSpectrum && (polymerPValue == (int)polymerPValue)) {
				spectraOutputStream << polymerPValue << ".0";
			} else {
				spectraOutputStream << polymerPValue;
			}

			// Base peak and second most intense
			spectraOutputStream << '\t'
				<< basePeakMz << (firstSpectrum && basePeakMz==0.0 ? ".0" : "") << '\t'
				<< basePeakIntensity << (firstSpectrum && basePeakIntensity==0.0 ? ".0" : "") <<'\t'
				<< secondPeakMz << (firstSpectrum && secondPeakMz==0.0 ? ".0" : "") <<'\t'
				<< secondPeakIntensity << (firstSpectrum && secondPeakIntensity==0.0 ? ".0" : "");

			// Status log
			spectraOutputStream << '\t'
				<< sourceCurrent << '\t'
				<< vacuumIonGauge << '\t'
				<< vacuumConvectronGauge << '\t'
				<< ftVacuumPenningGauge << '\t'
				<< ftVacuumPiraniGauge1 << '\t'
				<< ionMultiplier1 << '\t'
				<< ionMultiplier2 << '\t'
				<< ftCeMeasureVoltage << '\t'
				<< ftAnalyzerTemp;

			spectraOutputStream 
				<< '\n';
		}

		firstSpectrum=false;

		scan_num = fRawData->GetNextScanNum(scan_num);
	}

	if(!spectraFileName.empty()) {
		spectraOutputStream.close();
	}

	if(map!=NULL) {
		map->dumpEqualized(chromatogramMapFileName);
		delete map;
	}
}

void dumpToFile(std::string file, std::string contents) {
	std::ofstream out;
	out.exceptions(std::ofstream::failbit);
	out.open(file);
	out << std::setprecision(12);
	out << contents;
	out.close();
}

void extractTuneMethodData(Engine::Readers::FinniganRawData *fRawData, std::string fileName) {
	std::cout << "Extracting tune method information to file " << fileName << "." << std::endl;
	std::string data;
	fRawData->GetTuneMethod(&data);
	dumpToFile(fileName, data);	
}

void extractInstrumentMethodData(Engine::Readers::FinniganRawData *fRawData, std::string fileName) {
	std::cout << "Extracting instrument method information to file " << fileName << "." << std::endl;
	std::string data;
	fRawData->GetInstrumentMethod(&data);
	dumpToFile(fileName, data);	
}

void extractSampleInfoData(Engine::Readers::FinniganRawData *fRawData, std::string fileName) {
	std::cout << "Extracting sample information to file " << fileName << "." << std::endl;
	std::string data;
	fRawData->GetSampleInformation(&data);
	dumpToFile(fileName, data);	
}

void extractErrorLogData(Engine::Readers::FinniganRawData *fRawData, std::string fileName) {
	std::cout << "Extracting error log information to file " << fileName << "." << std::endl;
	std::string data;
	fRawData->GetErrorLog(&data);
	dumpToFile(fileName, data);	
}

int extractRawDataFile(std::string inputRawFileName, std::string infoFileName, std::string spectraFileName, std::string chromatogramMapFileName,
	std::string tuneMethodFileName, std::string instrumentMethodFileName, std::string sampleInfoFileName, std::string errorLogFileName) {
		using namespace std;

		clock_t startClock = clock();

		Engine::Readers::FinniganRawData * fRawData = NULL;

		cout << "Reading raw file: " << inputRawFileName << endl;

		try {
			fRawData = (Engine::Readers::FinniganRawData*)(Engine::Readers::ReaderFactory::GetRawData(Engine::Readers::FileType::FINNIGAN, const_cast<char*>(inputRawFileName.c_str())));		
		} catch(const char *exception) {
			cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " - " << exception << "\n";
			return 1;
		}

		try {
			if(!infoFileName.empty()) {
				extractInfoFile(fRawData, infoFileName.c_str());
			}

			if(!tuneMethodFileName.empty()) {
				extractTuneMethodData(fRawData, tuneMethodFileName);
			}

			if(!instrumentMethodFileName.empty()) {
				extractInstrumentMethodData(fRawData, instrumentMethodFileName);
			}

			if(!sampleInfoFileName.empty()) {
				extractSampleInfoData(fRawData, sampleInfoFileName);
			}

			if(!errorLogFileName.empty()) {
				extractErrorLogData(fRawData, errorLogFileName);
			}

			extractPerSpectrumData(fRawData, spectraFileName, chromatogramMapFileName);

		} catch(const char *exception) {
			std::cerr << "ERROR: problem generating output data file " << " - " << exception << " " << strerror(errno) << "\n";
			return 1;
		} catch(std::exception &e) {				
			std::cerr << "ERROR: problem generating output data file " << " - " << strerror(errno) << " - " << e.what() << "\n";
			return 1;
		}

		if(fRawData) {
			fRawData->Close();
			delete fRawData;
		}

		clock_t endClock = clock();

		std::cout << "Took " << ((double)endClock-startClock)/CLOCKS_PER_SEC << " seconds \n";

		return 0;
}

int extractRange(const char *inputRawFileName, const char *peaksFileName, double minMzValue, double maxMzValue) {
	clock_t startClock = clock();

	Engine::Readers::RawData * pRawData = NULL;

	try {
		pRawData = Engine::Readers::ReaderFactory::GetRawData(Engine::Readers::FileType::FINNIGAN, const_cast<char*>(inputRawFileName));		
	} catch(...) {
		std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << peaksFileName << "\n";
		return 1;
	}

	int errorCode=0;
	try {

		std::ofstream outputStream;
		outputStream.exceptions(std::ofstream::failbit);
		outputStream.open(peaksFileName);
		outputStream << std::setprecision(12);

		outputStream << ScanId << '\t' << RetentionTime << '\t' << Mz << '\t' << Intensity << '\n';

		std::vector<double> mzs;
		std::vector<double> intensities;

		for(int i=pRawData->GetFirstScanNum(); i<=pRawData->GetLastScanNum(); i=pRawData->GetNextScanNum(i)) {
			int msLevel = pRawData->GetMSLevel(i);
			// MS level 1 indicates a survey (FTMS) scan
			bool addThisScan = (msLevel==1);
			if(addThisScan) {			
				pRawData->GetRawData(mzs, intensities, i, minMzValue, maxMzValue);
				int scanNumber = i;
				double retentionTime = pRawData->GetScanTime(scanNumber);
				int numPeaks = mzs.size();
				for(int j=0; j<numPeaks; j++) {
					outputStream << scanNumber << '\t' << retentionTime << '\t' << mzs[j] << '\t' << intensities[j] << '\n';
				}
			}
		}		

		outputStream.close();
	} catch(...) {
		std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << peaksFileName << "\n";
		errorCode=1;
	}

	if(pRawData) {
		pRawData->Close();
		delete pRawData;
	}

	clock_t endClock = clock();

	std::cout << "Took " << ((double)endClock-startClock)/CLOCKS_PER_SEC << " seconds \n";

	return errorCode;

}

int generateMprcData (const char * inputRawFileName, const char * outputFileName, bool addMs2Peaks) {
	clock_t startClock = clock();

	Engine::Readers::RawData * pRawData = NULL;

	sqlite3 *db = NULL;
	char *zErrMsg = 0;
	int rc;

	try {
		pRawData = Engine::Readers::ReaderFactory::GetRawData(Engine::Readers::FileType::FINNIGAN, const_cast<char*>(inputRawFileName));		
	} catch(...) {
		std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << outputFileName << "\n";
		return 1;
	}

	int errorCode=0;
	try {
		// Open SQLite3 database
		rc = sqlite3_open(outputFileName, &db);
		if( rc ){
			fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
			sqlite3_close(db);
			exit(1);
		}

		pragmas(db);

		exec(db, "BEGIN;");
		createVersionTable(db);
		createContentTables(db);	
		if(errorCode==0) {
			errorCode = addSpectra(db, inputRawFileName, outputFileName, pRawData);
		}
		if(errorCode==0) {
			errorCode |= addMsScans(db, inputRawFileName, outputFileName, pRawData, addMs2Peaks);
		}
		exec(db, "COMMIT;");

	} catch(...) {
		std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << outputFileName << "\n";
		errorCode=1;
	}
	if(db!=NULL) {
		sqlite3_close(db);
	}

	if(pRawData) {
		pRawData->Close();
		delete pRawData;
	}

	clock_t endClock = clock();

	std::cout << "Took " << ((double)endClock-startClock)/CLOCKS_PER_SEC << " seconds \n";

	return errorCode;
}

std::string getOption(std::vector<std::string> & commandLineParams, std::string option) {

	using namespace std;

	vector<string>::iterator commandLineParamsItr = commandLineParams.begin();
	string value;

	while(commandLineParamsItr != commandLineParams.end()) {
		value = *commandLineParamsItr;
		if (value == option) {
			commandLineParamsItr++;

			if (commandLineParamsItr != commandLineParams.end()) {
				return *commandLineParamsItr;
			}
		}

		commandLineParamsItr++;
	}

	return "";
}

bool hasOption(std::vector<std::string> & commandLineParams, std::string option) {

	using namespace std;

	vector<string>::iterator commandLineParamsItr = commandLineParams.begin();
	string value;

	while(commandLineParamsItr != commandLineParams.end()) {
		value = *commandLineParamsItr;
		if (value == option) {
			return true;
		}

		commandLineParamsItr++;
	}

	return false;
}

std::string trim(std::string str) {
	using namespace std;

	int beginIndex = -1;
	int endIndex = -1;

	string trimedStr = str;

	const char* charStr = str.c_str();

	//Trim start of string
	for(unsigned int i = 0; i < trimedStr.length(); i++) {
		if (charStr[i] != ' ' && charStr[i] != '\t') {
			if (beginIndex == -1) {
				beginIndex = i;
			}

			endIndex = i + 1;
		}
	}

	if (beginIndex != -1) {
		trimedStr = str.substr(beginIndex, endIndex - beginIndex);
	}

	return trimedStr;
}

std::vector<std::string> getParamsFromFile(const char* inputParamsFile) {
	using namespace std;

	vector<string> paramVector;

	cout << "Reading parameter file " << inputParamsFile << "." << endl;

	ifstream inputStream;
	inputStream.open(inputParamsFile);

	if(inputStream.fail()) {
		cerr << "ERROR: problem reading parameter file " << " - " << strerror(errno) << "\n";
		return paramVector;
	}

	string line;

	while (!getline(inputStream, line).eof()) {
		paramVector.push_back(trim(line));
	}

	inputStream.close();

	return paramVector;
}

void printUsage() {
	std::cerr << "Usage: MprcExtractRaw <option> <parameters>\n";
	std::cerr << std::endl;
	std::cerr << "Options:" << std::endl;
	std::cerr << "  " << mprc << "\textract .mprc file" << std::endl;
	std::cerr << "  " << data << "\textract .tsv files for Swift QA" << std::endl;
	std::cerr << "  " << mzRange << "\textract .tsv files with peak data from given range" << std::endl;
	std::cerr << "  " << paramsFile << "\tobtain parameters from specified file"<< std::endl;
	std::cerr << std::endl;

	std::cerr << mprc << " parameters:" << std::endl;
	std::cerr << "  " << rawFile << " <thermo finnigan RAW file path> " << std::endl;
	std::cerr << "  " << outFile << " <output database (must not exist)>" << std::endl;
	std::cerr << "  " << ms2Switch << "  optional: enables extraction of peaks from ms2 spectra" << std::endl;
	std::cerr << std::endl;

	std::cerr << data << " parameters:" << std::endl;
	std::cerr << "  " << rawFile << " <thermo finnigan RAW file path> "  << std::endl;
	std::cerr << "  " << infoFile << " <info output file> " << std::endl;
	std::cerr << "\tRaw file specific information, such as instrument id, " << std::endl;
	std::cerr << "\toriginal raw file name, adquisition start time, and so on." << std::endl;
	std::cerr << "  " << spectraFile << " <spectra output file>" << std::endl;
	std::cerr << "\tSpectra information. Every row in the file represents a spectrum." << std::endl;
	std::cerr << "  " << chromatogramFile << " <chromatogram gif file>" << std::endl;
	std::cerr << "\tChromatogram as a gif image." << std::endl;
	std::cerr << "  " << tuneMethodFile << " <tune method file>" << std::endl;
	std::cerr << "\tTune method (copied verbatim as seen in XCalibur)" << std::endl;
	std::cerr << "  " << instrumentMethodFile << " <instrument method file>" << std::endl;
	std::cerr << "\tInstrument method (copied verbatim as seen in XCalibur)" << std::endl;
	std::cerr << "  " << sampleInformationFile << " <sample information file>" << std::endl;
	std::cerr << "\tSample information (copied verbatim as seen in XCalibur)" << std::endl;
	std::cerr << "  " << errorLogFile << " <error log file>" << std::endl;
	std::cerr << "\tError log. If there was no error, this file should be empty." << std::endl;

	std::cerr << std::endl;
	std::cerr << "Spectra file columns:" << std::endl;
	std::cerr << "* " << ScanId << "\t\tScan number" << std::endl;
	std::cerr << "* " << ParentMz << "\t\tThe m/z that was isolated to obtain MSn" << std::endl;
	std::cerr << "* " << TotalIonCurrent << "\t\t\tTotal Ion Current" << std::endl;
	std::cerr << "* " << RetentionTime << "\t\t\tRetention Time (minutes)" << std::endl;
	std::cerr << "* " << MsLevel << "\t\tMs Level (1-MS, 2-MS/MS, ...)" << std::endl;
	std::cerr << "* " << ParentScan << "\t\tID of the parent MS scan" << std::endl;
	std::cerr << "* " << ChildScans << "\t\tNumber of child scans (MS2) per MS scan" << std::endl;
	std::cerr << "* " << IonInjectionTimeMs << "\tIon injection time in ms" << std::endl;
	std::cerr << "* " << CycleTimeSeconds << "\t\tTime between two consecutive MS scans (seconds)" << std::endl;
	std::cerr << "* " << ElapsedTimeSeconds << "\t\tElapsed scan time (seconds)" << std::endl;
	std::cerr << "* " << DeadTimeSeconds << "\t\tDead time = Cycle time - Elapsed time (seconds)" << std::endl;
	std::cerr << "* " << TimeToNextScanSeconds << "\tRT(next scan)-RT(this scan) (seconds)" << std::endl;
	std::cerr << "* " << LockMassFound << "\t1 if lock mass was found in current scan" << std::endl;
	std::cerr << "* " << LockMassShift << "\tLock mass shift in PPM" << std::endl;
	std::cerr << "* " << DissociationType << "\tcid/etd" << std::endl;
	std::cerr << "* " << PolymerSegment << "\tSize of polymer segment (e.g. 44 Da for typical polymer)" << std::endl;
	std::cerr << "* " << PolymerOffset << "\tThe initial polymer mass (end of the polymer before segments start)" << std::endl;
	std::cerr << "* " << PolymerScore << "\tPolymer score" << std::endl;
	std::cerr << "* " << PolymerPValue << "\tProbability that a higher or equal polymer score could be achieved randomly" << std::endl;
	std::cerr << "* " << BasePeakMz << "\tBase peak m/z" << std::endl;
	std::cerr << "* " << BasePeakIntensity << "\tBase peak intensity" << std::endl;
	std::cerr << "* " << SecondPeakMz << "\tSecond highest peak m/z (>" << secondPeakMinDistanceFromBase << " Da from base peak)" << std::endl;
	std::cerr << "* " << SecondPeakIntensity << "\tSecond highest peak intensity" << std::endl;

	std::cerr << std::endl;

	std::cerr << mzRange << " parameters:" << std::endl;
	std::cerr << "  " << rawFile << " <thermo finnigan RAW file path> "  << std::endl;
	std::cerr << "  " << peaksFile << " <peaks data output file> " << std::endl;
	std::cerr << "  " << minMz << " <minimum M/Z> " << std::endl;
	std::cerr << "  " << maxMz << " <maximum M/Z> " << std::endl;
	std::cerr << std::endl;

	std::cerr << paramsFile << " <param file>" << std::endl;
	std::cerr << "\tThe params input file must contain all the command line " << std::endl;
	std::cerr << "\tparameters, one per line, including flags." << std::endl;
	std::cerr << "\tExample: " << std::endl;
	std::cerr << "\t\t" << data << std::endl;
	std::cerr << "\t\t" << rawFile << std::endl;
	std::cerr << "\t\t" << "<RAW file path>" << std::endl;
	std::cerr << "\t\t" << infoFile << std::endl;
	std::cerr << "\t\t" << "<info output file>" << std::endl;
	std::cerr << "\t\t" << spectraFile << std::endl;
	std::cerr << "\t\t" << "<spectra output file>" << std::endl;
}

void printOptionError(std::string option) {
	std::cerr << "Value for " << option << " must be specified." << std::endl;
	printUsage();
}

// MprcExtractRaw uses DeconMsn's libraries to extract mass/intensity pairs for survey spectra
// from Thermo Finnigan .RAW files (for now, we can support other formats like ABI as well).
// The program dumps all the extracted data into SqLite database in Mprc Data Format.
//
// See http://proteomics.mayo.edu/wiki/trac.cgi/wiki/MprcDataFormat for more information about the output format.
int main(int argc, char* argv[])
{
	if(argc < 3 || (argc >= 3 && (strcmp(argv[1], data.c_str())!=0 && strcmp(argv[1], mprc.c_str())!=0 && strcmp(argv[1], paramsFile.c_str())!=0))) {
		printUsage();
		return 1;
	}

	std::vector<std::string> paramVector;

	//Build command line parameter vector.
	if (strcmp(argv[1], paramsFile.c_str()) != 0) {
		for (int i = 0; i < argc; i++) {
			paramVector.push_back(argv[i]);
		}
	} else {
		paramVector = getParamsFromFile(argv[2]);
	}

	if (hasOption(paramVector, mprc)) {
		std::string inputRawFileName = getOption(paramVector, rawFile);
		if (inputRawFileName.empty()) {
			printOptionError(rawFile);
			return 1;
		}

		std::string outputFileName = getOption(paramVector, outFile);
		if (outputFileName.empty()) {
			printOptionError(outFile);
			return 1;
		}

		bool addMs2Peaks = hasOption(paramVector, ms2Switch);

		return generateMprcData (inputRawFileName.c_str(), outputFileName.c_str(), addMs2Peaks);

	} else if (hasOption(paramVector, data)) {
		std::string inputRawFileName = getOption(paramVector, rawFile);
		if (inputRawFileName.empty()) {
			printOptionError(rawFile);
			return 1;
		}

		std::string infoFileName = getOption(paramVector, infoFile);
		std::string spectraFileName = getOption(paramVector, spectraFile);
		std::string chromatogramMapFileName = getOption(paramVector, chromatogramFile);
		std::string tuneMethodFileName = getOption(paramVector, tuneMethodFile);
		std::string instrumentMethodFileName = getOption(paramVector, instrumentMethodFile);
		std::string sampleInfoFileName = getOption(paramVector, sampleInformationFile);
		std::string errorLogFileName = getOption(paramVector, errorLogFile);

		return extractRawDataFile(
			inputRawFileName, 
			infoFileName, 
			spectraFileName, 
			chromatogramMapFileName,
			tuneMethodFileName,
			instrumentMethodFileName,
			sampleInfoFileName,
			errorLogFileName);
	} else if (hasOption(paramVector, mzRange)) {
		std::string inputRawFileName = getOption(paramVector, rawFile);
		if (inputRawFileName.empty()) {
			printOptionError(rawFile);
			return 1;
		}

		std::string minMzStr = getOption(paramVector, minMz);
		if (minMzStr.empty()) {
			printOptionError(minMz);
			return 1;
		}
		double minMzValue = atof(minMzStr.c_str());

		std::string maxMzStr = getOption(paramVector, maxMz);
		if (maxMzStr.empty()) {
			printOptionError(maxMz);
			return 1;
		}
		double maxMzValue = atof(maxMzStr.c_str());

		std::string peaksFileName = getOption(paramVector, peaksFile);
		if (peaksFileName.empty()) {
			printOptionError(peaksFile);
			return 1;
		}

		return extractRange(inputRawFileName.c_str(), peaksFileName.c_str(), minMzValue, maxMzValue);
	}
	/*
	clock_t startClock = clock();

	Engine::Readers::RawData * pRawData = NULL;

	sqlite3 *db = NULL;
	char *zErrMsg = 0;
	int rc;

	try {
	pRawData = Engine::Readers::ReaderFactory::GetRawData(Engine::Readers::FileType::FINNIGAN, inputRawFileName);		
	} catch(...) {
	std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << outputFileName << "\n";
	return 1;
	}

	int errorCode=0;
	try {
	// Open SQLite3 database
	rc = sqlite3_open(outputFileName, &db);
	if( rc ){
	fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
	sqlite3_close(db);
	exit(1);
	}

	pragmas(db);

	exec(db, "BEGIN;");
	createVersionTable(db);
	createContentTables(db);	
	if(errorCode==0) {
	errorCode = addSpectra(db, inputRawFileName, outputFileName, pRawData);
	}
	if(errorCode==0) {
	errorCode |= addMsScans(db, inputRawFileName, outputFileName, pRawData, addMs2Peaks);
	}
	exec(db, "COMMIT;");

	} catch(...) {
	std::cerr << "ERROR: problem extracting raw data from " << inputRawFileName << " into " << outputFileName << "\n";
	errorCode=1;
	}
	if(db!=NULL) {
	sqlite3_close(db);
	}

	if(pRawData) {
	pRawData->Close();
	delete pRawData;
	}

	clock_t endClock = clock();

	std::cout << "Took " << ((double)endClock-startClock)/CLOCKS_PER_SEC << " seconds \n";

	return errorCode;
	*/
}

