MprcExtractRaw
==============

Extract data from Thermo .RAW files.

Utilizes modified NCRR's source code from http://ncrr.pnnl.gov/software/ to read the .RAW files.

The resulting data are saved either to:

* tab separated files
* .mprc SQLite database file
* .gif file for the chromatogram

Build
-----

We use Microsoft Visual Studio 2010. It should not be impossible to do the build using Visual C++ 2010 Express as well.

Usage
-----

    MprcExtractRaw <option> <parameters>

### Options
  --mprc        extract .mprc file

  --data        extract .tsv files for Swift QA

  --mzrange     extract .tsv files with peak data from given range

  --params      obtain parameters from specified file

#### --mprc parameters
  --raw <thermo finnigan RAW file path>

  --out <output database (must not exist)>

  --ms2  optional: enables extraction of peaks from ms2 spectra

#### --data parameters
  --raw <thermo finnigan RAW file path>

  --info <info output file>
        Raw file specific information, such as instrument id,
        original raw file name, adquisition start time, and so on.

  --spectra <spectra output file>
        Spectra information. Every row in the file represents a spectrum.

  --chromatogram <chromatogram gif file>
        Chromatogram as a gif image.

##### Spectra file columns

* Scan Id               Scan number
* Parent m/z            The m/z that was isolated to obtain MSn
* TIC                   Total Ion Current
* RT                    Retention Time (minutes)
* MS Level              Ms Level (1-MS, 2-MS/MS, ...)
* Parent Scan           ID of the parent MS scan
* Child Scans           Number of child scans (MS2) per MS scan
* Ion Injection Time    Ion injection time in ms
* Cycle Time            Time between two consecutive MS scans (seconds)
* Elapsed Time          Elapsed scan time (seconds)
* Dead Time             Dead time = Cycle time - Elapsed time (seconds)
* Time To Next Scan     RT(next scan)-RT(this scan) (seconds)
* Lock Mass Found       1 if lock mass was found in current scan
* Lock Mass Shift       Lock mass shift in PPM
* Dissociation Type     cid/etd
* Polymer Segment Size  Size of polymer segment (e.g. 44 Da for typical polymer)

* Polymer Offset        The initial polymer mass (end of the polymer before segments start)
* Polymer Score Polymer score
* Polymer p-value       Probability that a higher or equal polymer score could be achieved randomly

#### --mzrange parameters
  --raw <thermo finnigan RAW file path>
  --peaks <peaks data output file>
  --min <minimum M/Z>
  --max <maximum M/Z>

#### --params <param file>
        The params input file must contain all the command line
        parameters, one per line, including flags.
        Example:
                --data
                --raw
                <RAW file path>
                --info
                <info output file>
                --spectra
                <spectra output file>
