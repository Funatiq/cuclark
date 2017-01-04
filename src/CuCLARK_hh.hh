/*
   CuCLARK, CLARK for CUDA-enabled GPUs.
   Copyright 2016, Robin Kobus <rkobus@students.uni-mainz.de>
   
   based on CLARK version 1.1.3, CLAssifier based on Reduced K-mers.
   Copyright 2013-2016, Rachid Ounit <rouni001@cs.ucr.edu>


   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @author: Robin Kobus, masters student at Institute of Computer Science, JGU Mainz
 * @project: CuCLARK, Metagenomic Classification with CUDA-enabled GPUs
 * 
 * Changes:
 * Heavily modified database loading (CuClarkDB->read),
 * classification process (getObjectsDataComputeFullGPU)
 * and file output (printExtendedResultsSynced).
 * Kept database creation process of CLARK with exception of using k-mers in canonical form.
 * Many smaller changes.
 */

#ifndef CUCLARK_HH
#define CUCLARK_HH

#include<cstring>
#include<vector>
#include<iostream>
#include<sstream>
#include<cstdlib>
#include "./dataType.hh"
#include "./HashTableStorage_hh.hh"
#include "./CuClarkDB.cuh"

//~ #include <bitset>

#define MAXRSIZE	10000

template <typename HKMERr>
class CuCLARK
{
	private:
		std::vector< std::pair< std::string, std::string > >	m_targetsID;
		std::vector< std::string >				m_labels;
		std::vector< std::string >				m_labels_c;
		std::vector< std::string >				m_targetsName;

		// Using fasta/fastq files
		ITYPE						m_nbObjects;

		bool						m_isFastaFile;
		
		// batch data
		std::vector< std::vector<size_t> >		m_seqSNames;	// name starting position
		std::vector< std::vector<size_t> >      m_seqENames;	// name starting position
		std::vector< std::vector<size_t> >    	m_readsSPos; 	// read starting position
		std::vector< std::vector<size_t> >    	m_readsEPos; 	// read ending position
		std::vector< std::vector<size_t> >    	m_readsLength;
		std::vector< size_t >					m_posReads;
//// new	
		vector<bool> 							m_batchScheduled;

		// Options for loading db				
		bool					m_isLightLoading;
//// new
		size_t					m_dbParts;

		const char*				m_folder;
		bool					m_isSpecificTargetsPresent;

		// Dictionary k-mers <---> TargetsID
		EHashtable<HKMERr, rElement> *	m_centralHt;
//// new db		
		CuClarkDB<HKMERr> * 			m_cuClarkDb;

//// new
		// Tables for storing results
		RESULTS*		m_fullResults;
		RESULTS*		m_finalResults;
		size_t			m_resultRowSize;
		size_t			m_finalResultsRowSize;
		
		size_t			m_nbCPU;
//// new		
		size_t			m_numDevices;
		size_t			m_numBatches;
		
		const size_t	m_kmerSize;
		const uint8_t	m_k;
		const ITYPE		m_minCountTarget;
		uint64_t		m_iterKmers;
		ITYPE			m_minCountObject;
		bool			m_isExtended;
		bool			m_isPaired;

		// Tables storing common and repetitive values for reading sequences
		int					m_Letter[256];
		int					m_table[256];
		int					m_rTable[256];
		uint64_t		    		m_pTable[4];
		std::vector< std::vector<uint64_t> >	m_powerTable;
		int	 				m_separators[256];

	public:
		CuCLARK(const size_t& 		_kmerLength,
				const char* 		_filesName,
				const char* 		_folderName,
				const ITYPE& 		_minCountT 	= 0,
				const bool&			_creatingkmfiles= false,
				const bool&			_isLightLoading = false,
				const uint64_t&		_iterKmers	= 2,
				const size_t& 		_nbCPU		= 1,
				const ITYPE&		_samplingFactor = 1,
				const bool& 		_mmapLoading 	= false,
				const size_t&		_numBatches = 1,
				const size_t&		_numDevices = 1
		     );

		~CuCLARK();

		void runSimple(const char* 		_fileTofilesname, 
				const char* 		_fileResult,
				const ITYPE& 		_minCountO 	= 0
			      );

		void run(const char*	_filesToObjects,
				const char* 	_fileToResults,
				const ITYPE& 	_minCountO 	= 0,
				const bool&		_isExtended	= false
			);

		void run(const char*	_pairedfile1,
				const char*		_pairedfile2,
				const char*		_fileToResults,
				const ITYPE&	_minCountO  = 0,
				const bool&		_isExtended = false
			);

		void clear();

	private:

		void createTargetFilesNames(std::vector< std::string >& 	_filesHT, 
				std::vector< std::string >& 			_filesHTC
				) const;

		void loadSpecificTargetSets(const std::vector<std::string>& 	_filesHT, 
				const std::vector<std::string>& 		_filesHTC, 
				const size_t& 					_sizeMotherHT,
				const ITYPE&            			_samplingFactor, 
				const bool&  					_mmapLoading
				);

		size_t makeSpecificTargetSets(const std::vector<std::string>& 	_filesHT, 
				const std::vector<std::string>& 		_filesHTC
				) const;

		void getObjectsDataComputeFullGPU(const uint8_t *			_map,
				const size_t& 					nb,
				const char* _fileResult
				);
				
		bool getTargetsData(const char* 				_filesName, 
				std::vector< std::string >&		 	_filesHT, 
				std::vector< std::string >&		 	_filesHTC,
				const bool&					_creatingkmfiles,
				const ITYPE&					_samplingfactor = 1
				);

		void print(const bool& 						_creatingkmfiles = false,
				const ITYPE&                                    _samplingfactor = 1
			  ) const;
				
		void printExtendedResultsSynced(const uint8_t * _map,  const char* _fileResult);
		
		void printSpeedStats(const struct timeval& 			_requestEnd, 
				const struct timeval& 				_requestStart, 
				const char* 					_fileResult
				)  const;

		void getdbName(char * 						_dbname
			      ) const;
};

#endif


#ifdef _OPENMP
#include <omp.h>
#endif

#include <string.h>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "./file.hh"
#include "./kmersConversion.hh"
#include "./analyser.hh"

using namespace std;

template <typename HKMERr>
CuCLARK<HKMERr>::CuCLARK(const size_t& 	_kmerLength,
		const char* 		_filesName,
		const char* 		_folderName,
		const ITYPE& 		_minCountT,
		const bool&      	_creatingkmfiles,
		const bool&     	_isLightLoading,
		const uint64_t&		_iterKmers,
		const size_t& 		_nbCPU,
		const ITYPE&        _samplingFactor,	
		const bool&     	_mmapLoading,
		const size_t&		_numBatches,
		const size_t&		_numDevices
		): 
	m_nbCPU(_nbCPU), 
	m_kmerSize(_kmerLength), m_k((uint8_t) _kmerLength),
	m_nbObjects(0),
	m_folder(_folderName), 
	m_minCountTarget(_minCountT), 
	m_powerTable(m_kmerSize),
	m_isFastaFile(true),
	m_iterKmers(_iterKmers),
	m_isLightLoading(_isLightLoading),
	m_posReads(_nbCPU),
	m_isPaired(false),
	m_numBatches(_numBatches),
	m_numDevices(_numDevices)
{

#ifdef _OPENMP
	omp_set_num_threads(m_nbCPU);
#else
	m_nbCPU = 1;
#endif
	m_posReads.resize(m_numBatches);
	m_readsLength.resize(m_numBatches); 
	m_readsEPos.resize(m_numBatches);
	m_readsSPos.resize(m_numBatches);
	m_seqENames.resize(m_numBatches);
	m_seqSNames.resize(m_numBatches);
	
	m_batchScheduled.resize(m_numBatches);

	size_t base = 1;
	for(size_t p = 0; p < m_kmerSize ; p++)
	{       
		m_powerTable[p].resize(4);
		m_powerTable[p][0] = 0; m_powerTable[p][1] = base; m_powerTable[p][2] = base*2; m_powerTable[p][3] = base*3;  
		base *= 4;
	}
	m_pTable[0] = m_powerTable[m_kmerSize-1][0]; 
	m_pTable[1] = m_powerTable[m_kmerSize-1][1];
	m_pTable[2] = m_powerTable[m_kmerSize-1][2];
	m_pTable[3] = m_powerTable[m_kmerSize-1][3];

	for(size_t t = 0; t < 256 ; t++)
	{	m_table[t] = -1; m_separators[t] = 0; m_rTable[t] = -1;	m_Letter[t] = -1;}

	// Letter definition (Fastq files)
	for(size_t t = 65; t < 91; t++)
	{	m_Letter[t] = 0;}
	for(size_t t = 97; t < 123; t++)
	{       m_Letter[t] = 0;}
	// Nucleotide definitions (DNA)
	m_table['A'] = 0;	m_table['C'] = 1;       m_table['G'] = 2;	m_table['T'] = 3;
	m_table['a'] = 0;       m_table['c'] = 1;       m_table['g'] = 2;       m_table['t'] = 3;
	// Nucleotide definitions (RNA)
	m_table['U'] = 3; 	m_table['u'] = 3;

	// Other symbol definitions
	m_table['>'] = -2;	m_table['\n']= -10;
	// Reverse - Nucleotide definitions (DNA)
	m_rTable['A'] = 3;      m_rTable['C'] = 2;       m_rTable['G'] = 1;       m_rTable['T'] = 0;
	m_rTable['a'] = 3;      m_rTable['c'] = 2;       m_rTable['g'] = 1;       m_rTable['t'] = 0;
	// Reverse - Nucleotide definitions (RNA)
	m_rTable['U'] = 0;      m_rTable['u'] = 0;

	// Reverse - Other symbol definitions
	m_rTable['>'] = -2;     m_rTable['\n']= -10;

	m_separators[' '] = 1;	m_separators['\t'] = 1; 	m_separators['\n'] = 1;

	vector<string> filesHT, filesHTC;
	size_t sizeMotherHT = 0;
	if (!getTargetsData(_filesName, filesHT, filesHTC, _creatingkmfiles, _samplingFactor))
	{
		cerr << "Starting the creation of the database of targets specific " << m_kmerSize << "-mers from input files..." << endl;
		sizeMotherHT = makeSpecificTargetSets(filesHT, filesHTC);
	}
	loadSpecificTargetSets(filesHT, filesHTC, sizeMotherHT, _samplingFactor, _mmapLoading);
}

template <typename HKMERr>
CuCLARK<HKMERr>::~CuCLARK()
{
	clear();
	delete m_centralHt;
}

template <typename HKMERr>
void CuCLARK<HKMERr>::clear()
{
	m_nbObjects = 0;

	//~ m_batchScheduled.clear();
	for(size_t  i = 0; i < m_numBatches; i++)
	{
		m_seqSNames[i].clear();
		m_seqENames[i].clear();
		m_readsEPos[i].clear();
		m_readsSPos[i].clear();
		m_readsLength[i].clear();
		
		m_batchScheduled[i] = false;
	}
	
	m_cuClarkDb->free();
}

/**
 *  Create targets specific k-mers files.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::createTargetFilesNames(vector<string>& _filesHT, vector<string>& _filesHTC) const
{
	for(size_t t = 0 ; t < m_labels_c.size(); t++)
	{
		char * fname = (char*) calloc(100, sizeof(char));
		if (m_isLightLoading)
		{
			sprintf(fname, "%s/%s_k%lu_light.ht", m_folder, m_labels_c[t].c_str(), (size_t) m_kmerSize);
		}
		else
		{
			sprintf(fname, "%s/%s_k%lu.ht", m_folder, m_labels_c[t].c_str(),(size_t) m_kmerSize);
		}
		string nameHTO(fname);
		_filesHTC.push_back(nameHTO);
		free(fname);
		fname = NULL;
	}

	for(size_t t = 0 ; t < m_labels.size(); t++)
	{
		char * fname = (char*) calloc(100, sizeof(char));
		if (m_isLightLoading)
		{
			sprintf(fname, "%s/%s_k%lu_light.ht", m_folder, m_labels_c[t].c_str(), (size_t) m_kmerSize);
		}
		else
		{
			sprintf(fname, "%s/%s_k%lu.ht", m_folder, m_labels[t].c_str(), (size_t) m_kmerSize);
		}
		string nameHT(fname);
		_filesHT.push_back(nameHT);
		free(fname);
		fname = NULL;
	}
}

/**
 * Check files and set up to run classification.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::run(const char* _filesToObjects, const char* _fileToResults, const ITYPE& _minCountO, const bool& _isExtended)
{
	FILE* fd = fopen(_fileToResults, "r");
	m_isPaired = false;
	m_isExtended = _isExtended;
	if (fd == NULL )
	{
		cout << "Processing file \'" << _filesToObjects << "\' in " << m_numBatches << " batches using "<< m_nbCPU << " CPU thread(s)." <<  endl;
		CuCLARK::runSimple(_filesToObjects, _fileToResults, _minCountO);
		return;
	}
	fclose(fd);

	fd = fopen(_filesToObjects, "r");
	string line = "";
	getLineFromFile(fd, line);
	vector<string> ele;
	vector<char> seps;
	seps.push_back(' ');
	seps.push_back('\t');
	seps.push_back(',');
	fclose(fd);

	getElementsFromLine(line, seps, ele);	

	if (line[0] == '>' || line[0] == '@' || ele.size() == 2)
	{
		cout << "Processing file\'" << _filesToObjects << "\' in " << m_numBatches << " batches using "<< m_nbCPU << " CPU thread(s)." <<  endl;
		CuCLARK::runSimple(_filesToObjects, _fileToResults, _minCountO);
		return;
	}
	
	// run for multiple inputs
	FILE * r_fd = fopen(_fileToResults, "r");
	FILE * o_fd = fopen(_filesToObjects, "r");
	string o_line = "", r_line = "";
	cout <<  "Using " << omp_get_max_threads() << " CPU thread(s)." <<  endl;
	while (getLineFromFile(o_fd, o_line) && getLineFromFile(r_fd, r_line))
	{
		cout << "> Processing file \'" << o_line.c_str() << "\' in " << m_numBatches << " batches."<< endl;
		CuCLARK::runSimple(o_line.c_str(), r_line.c_str(), _minCountO); 
	}
	fclose(r_fd); 
	fclose(o_fd);
	return;
}

/**
 * Check files and set up to run classification (paired version).
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::run(const char* _pairedfile1, const char* _pairedfile2, const char* _fileToResults, const ITYPE& _minCountO, const bool& _isExtended)
{
	FILE* fd 	= fopen(_fileToResults, "r");
	m_isPaired 	= true;
	m_isExtended 	= _isExtended;
	char * mergedFiles = NULL;
	if (fd == NULL )
	{
		// Merge _pairedfile1 + _pairedfile2
		mergedFiles = (char *) calloc(strlen(_pairedfile1)+25, 1);
		sprintf(mergedFiles,"%s_ConcatenatedByCLARK.fa",_pairedfile1);
		mergePairedFiles(_pairedfile1, _pairedfile2, mergedFiles);
		cout << "Processing file: \'" << mergedFiles << "\' in " << m_numBatches << " batches using "<< m_nbCPU << " CPU thread(s)." <<  endl;
		CuCLARK::runSimple(mergedFiles, _fileToResults, _minCountO);
		// Delete file
		deleteFile(mergedFiles);
		free(mergedFiles);
		mergedFiles = NULL;
		return;
	}
	fclose(fd);

	fd = fopen(_pairedfile1, "r");
	string line 	= "";
	getLineFromFile(fd, line);
	vector<string> ele;
	vector<char> seps;
	seps.push_back(' ');
	seps.push_back('\t');
	seps.push_back(',');
	fclose(fd);

	getElementsFromLine(line, seps, ele);
	if (line[0] == '>' || line[0] == '@' || ele.size() == 2)
	{
		// Merge _pairedfile1 + _pairedfile2 
		mergedFiles = (char *) calloc(strlen(_pairedfile1)+25, 1);
		sprintf(mergedFiles,"%s_ConcatenatedByCLARK.fa",_pairedfile1);
		mergePairedFiles(_pairedfile1, _pairedfile2, mergedFiles);
		cout << "Processing file: \'" << mergedFiles << "\' in " << m_numBatches << " batches using "<< m_nbCPU << " CPU thread(s)." <<  endl;
		CuCLARK::runSimple(mergedFiles, _fileToResults, _minCountO);
		// Delete file
		deleteFile(mergedFiles);
		free(mergedFiles);
		mergedFiles = NULL;
		return;
	}
	
	// run for multiple inputs
	FILE * r_fd 	= fopen(_fileToResults, "r");
	FILE * o1_fd 	= fopen(_pairedfile1, "r");
	FILE * o2_fd 	= fopen(_pairedfile2, "r");
	string o1_line 	= "", o2_line   = "", r_line = "";
	cout <<  "Using " << omp_get_max_threads() << " CPU thread(s)." <<  endl;
	while (getLineFromFile(o1_fd, o1_line) && getLineFromFile(o2_fd, o2_line) && getLineFromFile(r_fd, r_line))
	{
		// Merge _pairedfile1 + _pairedfile2 
		mergedFiles = (char *) calloc(strlen(o1_line.c_str())+25, 1);
		sprintf(mergedFiles,"%s_ConcatenatedByCLARK.fa",o1_line.c_str());
		mergePairedFiles(o1_line.c_str(), o2_line.c_str(), mergedFiles);

		cout << "> Processing file: \'" << mergedFiles << "\' in " << m_numBatches << " batches."<< endl;
		CuCLARK::runSimple(mergedFiles, r_line.c_str(), _minCountO);
		// Delete file
		deleteFile(mergedFiles);
		free(mergedFiles);
		mergedFiles = NULL;
	}
	fclose(r_fd);
	fclose(o1_fd);
	fclose(o2_fd);
	return;
}

/**
 * Run the classification for a simgle input file.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::runSimple(const char* _fileTofilesname, const char* _fileResult, const ITYPE& _minCountO)
{
	clear();
	m_cuClarkDb->swapDbParts();
	m_cuClarkDb->sync();
	
	m_isFastaFile		= true;
	m_minCountObject 	= _minCountO;

	size_t fileSize;
	std::ifstream in(_fileTofilesname, std::ios::binary | std::ios::ate);
	fileSize = in.tellg();

	int fd = open(_fileTofilesname, O_RDONLY);
	if (fd == -1 || fileSize == 0)
	{
		cerr << "Failed to open " << _fileTofilesname << endl;
		return;
	}
	uint8_t *map;
	map = (uint8_t*) mmap(0, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
	if ( map == MAP_FAILED )
	{
		close(fd);
		cerr << "Failed to mmapping the file." << endl;
		return;
	}
	// Checking file to store result:
	string sfileResult(_fileResult);
	sfileResult += ".csv";
	const char* fileResult = sfileResult.c_str();
	// Try to access and erase content of the file If non-empty
	FILE * _fout = fopen(fileResult,"w");
	if (_fout == NULL)
	{
		cerr << "Failed to create/open file result: " << fileResult << endl;
		return;
	}

	struct timeval requestStart, requestEnd;

	fclose(_fout);

	if (true)
	{
		gettimeofday(&requestStart, NULL);	
		///////////////////////////////////////////////////////////////////////
		//~ getObjectsDataComputeFull(map, fileSize);
		getObjectsDataComputeFullGPU(map, fileSize, fileResult);
		//~ printExtendedResults(fileResult);
		///////////////////////////////////////////////////////////////////////
		gettimeofday(&requestEnd, NULL);
		// Measurement execution time
		printSpeedStats(requestEnd,requestStart,fileResult);

		msync(map, fileSize, MS_SYNC);
		if (munmap(map, fileSize) == -1)
		{       cerr << "Error un-mmapping the file." << endl;}
		close(fd);
		return;
	}

	msync(map, fileSize, MS_SYNC);
	if (munmap(map, fileSize) == -1)
	{       cerr << "Error un-mmapping the file." << endl;;}
	close(fd);

	return;
}

/**
 *  Get database name composed of database parameters.
 * 	This enables to store different databases in the same folder.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::getdbName(char *   _dbname) const
{
	size_t sizeHTS =  m_labels.size() + m_labels_c.size(); 
	if (m_isLightLoading)
	{
		sprintf(_dbname,"%s/db_central_k%lu_t%lu_s%lu_m%lu_light_%lu.tsk",m_folder,(size_t)m_kmerSize,sizeHTS,(size_t) HTSIZE,(size_t)m_minCountTarget,(size_t) m_iterKmers);	
	}
	else
	{
		sprintf(_dbname,"%s/db_central_k%lu_t%lu_s%lu_m%lu.tsk",m_folder,(size_t)m_kmerSize,sizeHTS,(size_t) HTSIZE,(size_t)m_minCountTarget);
	}
}

/**
 * Load database of target-specific k-mers.
 * If that fails, try to recover it from saved targets-specific data.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::loadSpecificTargetSets(const vector<string>& _filesHT, 
		const vector<string>& 	_filesHTC, 
		const size_t& 		_sizeMotherHT,
		const ITYPE& 		_samplingFactor,	
		const bool&  		_mmapLoading
		)
{
	size_t kmersLoaded = 0;
	ITYPE minCount = m_minCountTarget;
////	m_centralHt = new EHashtable<HKMERr, rElement>(m_kmerSize, m_labels, m_labels_c);
	m_cuClarkDb = new CuClarkDB<HKMERr>(m_numDevices, m_kmerSize, m_numBatches, m_targetsName.size()-1);
	char * cfname = (char*) calloc(130, sizeof(char));
	getdbName(cfname);

	cerr << "Loading database [" << cfname << ".*] (s=" << _samplingFactor<< ")..." << endl;

	size_t fileSize;
	
#ifdef TIME_DBLOADING
	struct timeval requestStart, requestEnd;
	gettimeofday(&requestStart, NULL);
#endif
////	if (m_centralHt->Read(cfname, fileSize, m_nbCPU, _samplingFactor, _mmapLoading))
	if (m_cuClarkDb->read(cfname, fileSize, m_dbParts, _samplingFactor, _mmapLoading))
	{
#ifdef TIME_DBLOADING
		gettimeofday(&requestEnd, NULL);
		double diff = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_usec - requestStart.tv_usec) / 1000000.0;
		cerr << "Loading time: " << diff << " s\n";
#endif
		//~ cerr << "Loading done (database size: " << fileSize / 1000000<<" MB)" << endl;
		free(cfname);
		cfname = NULL;
		return;
	}
	if ( _filesHT.size() + _filesHTC.size() == 0)
	{
		cerr << "Failed to find the database." << endl;
		exit(-1);	
	}
	cerr << "The database will be recovered from saved targets-specific data." << endl;

	m_centralHt = new EHashtable<HKMERr, rElement>(m_kmerSize, m_labels, m_labels_c);
	
	for(size_t t = 0 ; t < _filesHT.size(); t++)
	{
		string nameHT =  _filesHT[t];

		FILE* fd = fopen(nameHT.c_str(),"r");
		if (fd == NULL)
		{
			cerr << "Failed to open " << nameHT << endl;
		}
		else
		{
			fclose(fd);
			m_centralHt->Load(nameHT, m_labels[t], minCount);
			kmersLoaded = m_centralHt->Size();
		}
		cerr << "\rDataset " << t+1 << " loaded.   " ;
	}
	for(size_t t = 0 ; t < _filesHTC.size(); t++)
	{
		string nameHTO =  _filesHTC[t];

		FILE* fd = fopen(nameHTO.c_str(),"r");
		if (fd == NULL)
		{
			cerr << "Failed to open " << nameHTO << endl;
		}
		else
		{
			fclose(fd);
			m_centralHt->Load(nameHTO, m_labels_c[t], minCount);
			kmersLoaded = m_centralHt->Size();
		}
		cerr << "\rDataset " << t + 1 + _filesHT.size() << " loaded.    " ;
	}
	cerr << kmersLoaded  << " " << m_kmerSize << "-mers finally loaded. ";
	cerr << "Creating database in disk..." << endl;

	m_centralHt->SortAllHashTable();
	m_centralHt->Write(cfname, 0, false);
	free(cfname);
	cfname = NULL;
	cerr << "Central Hashtable successfully stored in disk." << endl;
	exit(-1);
}

/**
 *  Create database of targets specific k-mers.
 */
template <typename HKMERr>
size_t CuCLARK<HKMERr>::makeSpecificTargetSets(const vector<string>& _filesHT, const vector<string>& _filesHTC) const
{
	size_t nt = 0;
	if (m_isLightLoading)
	{
		EHashtable<HKMERr, lElement> commonKmersHT(m_kmerSize, m_labels, m_labels_c);
		// read kmers from files and add to HT
		for(size_t t = 0 ; t < m_targetsID.size(); t++)
		{
			FILE* fd = fopen(m_targetsID[t].first.c_str(),"r");
			if (fd == NULL)
			{        cerr << "Failed to open " << m_targetsID[t].first << endl;
				continue;
			}
			char c[MAXRSIZE];
			size_t len = fread(&c, 1, MAXRSIZE, fd), i = 0;
			uint64_t _km_f = 0, _km_r = 0;		// kmer, reverse kmer
			uint8_t cpt = 0;					// kmer length
			uint64_t iter = 0;
			// fasta
			if (len > 0 && c[0] == '>')
			{
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)		// process read
						{
							nt++;
							_km_r <<= 2;    
							_km_r += 3 - m_table[c[i]];
							if (cpt  == m_kmerSize - 1)		// kmer complete
							{
								//[/////////////////////////////////////////////////////////////////////
								if ( iter % m_iterKmers == 0  )		// skip kmers for light loading
								{       commonKmersHT.addElement(_km_r, m_targetsID[t].second, 1);}
								_km_r = 0; cpt = 0; 	// reset kmer
								i++;
								iter++;
								continue;
								//]//////////////////////////////////////////////////////////////////////
							}
							cpt++;  
							i++;
							continue;
						}
						if (m_table[c[i]] == -10)	// c[i] == '\n', skip
						{       i++;
							continue;   }
						if (m_table[c[i]] == -1)	// c[i] == unknown nucleotide, reset kmer, skip
						{        nt++; _km_r = 0; cpt = 0; i++;
							continue;  }
						if (m_table[c[i]] == -2)	// c[i] == '>', skip to next line
						{
							_km_r = 0; cpt = 0; 
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")    ";
				continue;
			}
			// fastq
			if (len > 0 && c[0] == '@')
			{
				while (len != 0)
				{
					while(i < len && c[i] != '\n')
					{       i++;}
					if (i < len)
					{   break;}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				i++;
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							_km_r <<= 2; 
							_km_r += 3 - m_table[c[i]];
							if (cpt  == m_kmerSize - 1)
							{
								//[/////////////////////////////////////////////////////////////////////
								if ( iter % m_iterKmers == 0 )
								{       commonKmersHT.addElement(_km_r, m_targetsID[t].second, 1);}
								_km_r = 0; cpt = 0; 
								i++;
								iter++;
								continue;
								//]//////////////////////////////////////////////////////////////////////
							}
							cpt++;
							i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{
							_km_r = 0; cpt = 0;
							// Pass third line
							i++;
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass fouth line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass first line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						if (m_table[c[i]] == -1)
						{
							nt++;
							_km_r = 0; cpt = 0; 
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
				continue;
			}
			// spectrum form
			fseek(fd, 0, 0);
			string s_kmer = "";
			ITYPE val;
			uint8_t counter = 0;
			while (getFirstAndSecondElementInLine(fd, s_kmer, val))
			{
				if (counter % m_iterKmers == 0 && val > m_minCountTarget)
				{
					commonKmersHT.addElement(s_kmer, m_targetsID[t].second, (size_t) val);
					counter = 0;
				}
				counter++;
				continue;
			}
			fclose(fd);
			cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
		}
		cerr << nt << " nt read in total." << endl;
		cerr << "Mother Hashtable successfully built. "<<commonKmersHT.Size()<<" " << m_kmerSize << "-mers stored." <<  endl;
		size_t sizeMotherTable = commonKmersHT.Size();

		commonKmersHT.SortAllHashTable(2);
		commonKmersHT.RemoveCommon(m_labels_c, m_minCountTarget);
		char * cfname = (char*) calloc(130, sizeof(char));
		getdbName(cfname);

		cerr << "Creating light database in disk..." << endl;
		uint64_t nbElement = commonKmersHT.Write(cfname,2);
		free(cfname);
		cfname = NULL;
		cerr << nbElement << " " << m_kmerSize << "-mers successfully stored in database." << endl;

		return sizeMotherTable;
	}
	if (_filesHTC.size() + _filesHT.size() == 0)
	{
		EHashtable<HKMERr, lElement> commonKmersHT(m_kmerSize, m_labels, m_labels_c);
		for(size_t t = 0 ; t < m_targetsID.size(); t++)
		{
			FILE* fd = fopen(m_targetsID[t].first.c_str(),"r");
			if (fd == NULL)
			{
				cerr << "Failed to open " << m_targetsID[t].first << endl;
			}
			else
			{
				char c[MAXRSIZE];
				size_t len = fread(&c, 1, MAXRSIZE, fd), i = 0;
				uint64_t _km_f = 0, _km_r = 0;
				bool _isfull = false;
				uint8_t cpt = 0;		
				// fasta
				if (len > 0 && c[0] == '>')
				{
					while (len != 0)
					{
						while( i < len)
						{
							if (m_table[c[i]] >= 0)
							{
								nt++;
								if (_isfull)
								{
									_km_f >>= 2;
									_km_f += m_powerTable[cpt][m_table[c[i]]];

									/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_f, m_targetsID[t].second, 1);
									///////////////////////////////////////////////////////////////////////
									i++;
									continue;
								}
								_km_r <<= 2;    
								_km_r += 3 - m_table[c[i]];
								if (cpt  == m_kmerSize -1)
								{       
									_isfull = true;
									//[/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_r, m_targetsID[t].second, 1);
									//]//////////////////////////////////////////////////////////////////////
									_km_f = _km_r;
									// The following 6 lines come from Jellyfish source code
									_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
									_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
									_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
									_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
									_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
									_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
									i++;
									continue;
								}
								cpt++;  i++;    
								continue;
							}
							if (m_table[c[i]] == -10)
							{	i++;
								continue;   }
							if (m_table[c[i]] == -1)
							{	nt++; _km_r = 0; cpt = 0; _isfull = false; i++; 
								continue;  }
							if (m_table[c[i]] == -2)
							{
								_km_r = 0; cpt = 0; _isfull = false; 
								while (len != 0)
								{
									while(i < len && c[i] != '\n')  {	i++;}
									if (i < len)		{   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;   
								continue;
							}
							cerr << m_targetsID[t].first << ": " ;
							cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
						}
						len = fread(&c, 1, MAXRSIZE, fd );
						i = 0;
					}
					fclose(fd);
					cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")    ";
					continue;	
				}
				// fastq
				if (len > 0 && c[0] == '@')
				{
					i = 0;
					while (len != 0)
					{
						while(i < len && c[i] != '\n')            {	i++;}
						if (i < len)      {   break;}
						len = fread(&c, 1, MAXRSIZE, fd );
						i = 0;
					}
					i++;
					while (len != 0)
					{
						while( i < len)
						{
							if (m_table[c[i]] >= 0)
							{
								nt++;
								if (_isfull)
								{
									_km_f >>= 2; 
									_km_f += m_powerTable[cpt][m_table[c[i]]];
									/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_f, m_targetsID[t].second, 1);
									///////////////////////////////////////////////////////////////////////
									i++;    
									continue;
								}
								_km_r <<= 2;    
								_km_r += 3 - m_table[c[i]];
								if ( cpt  == m_kmerSize - 1)
								{       _isfull = true;
									//[/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_r, m_targetsID[t].second, 1);
									//]//////////////////////////////////////////////////////////////////////
									_km_f = _km_r;
									// The following 6 lines come from Jellyfish source code
									_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
									_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
									_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
									_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
									_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
									_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
									i++;
									continue;
								}
								cpt++;  
								i++;    
								continue;
							}
							if (m_table[c[i]] == -10)
							{
								_km_r = 0; cpt = 0; _isfull = false; 
								// Pass third line
								i++;
								while (len != 0)
								{
									while(i < len && c[i] != '\n')		{	i++;}
									if (i < len)   		{   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;
								// Pass fouth line
								while (len != 0)
								{
									while(i < len && c[i] != '\n')         {	i++;}
									if (i < len)         {   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;
								// Pass first line
								while (len != 0)
								{
									while(i < len && c[i] != '\n')        {		i++;}
									if (i < len){   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;    
								continue;
							}
							if (m_table[c[i]] == -1)
							{ 	nt++;   
								_km_r = 0; cpt = 0; _isfull = false;
								i++;    
								continue;
							}
							cerr << m_targetsID[t].first << ": " ;
							cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
						}
						len = fread(&c, 1, MAXRSIZE, fd );
						i = 0;  
					}
					fclose(fd);
					cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
					continue;
				}
				// spectrum form
				fseek(fd, 0, 0);
				string s_kmer = "";
				ITYPE val;
				while (getFirstAndSecondElementInLine(fd, s_kmer, val))
				{	if (val > m_minCountTarget)
					{	commonKmersHT.addElement(s_kmer, m_targetsID[t].second, (size_t) val);}
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
			}
		}
		cerr << nt << " nt read in total." << endl;
		cerr << "Mother Hashtable successfully built. " << commonKmersHT.Size() << " " << m_kmerSize << "-mers stored." <<  endl;
		size_t sizeMotherTable = commonKmersHT.Size();

		commonKmersHT.SortAllHashTable(2);
		commonKmersHT.RemoveCommon(m_labels_c, m_minCountTarget);
		char * cfname = (char*) calloc(130, sizeof(char));
		getdbName(cfname);
		cerr << "Creating database in disk..." << endl;	
		uint64_t nbElement = commonKmersHT.Write(cfname,2);
		free(cfname);
		cfname = NULL;
		cerr << nbElement << " " << m_kmerSize << "-mers successfully stored in database." << endl;

		return sizeMotherTable;
	}
	///////////////////////////////////////////////////////////////////////////////
	EHashtable<HKMERr, Element> commonKmersHT(m_kmerSize, m_labels, m_labels_c);
	for(size_t t = 0 ; t < m_targetsID.size(); t++)
	{
		FILE* fd = fopen(m_targetsID[t].first.c_str(),"r");
		if (fd == NULL)
		{
			cerr << "Failed to open " << m_targetsID[t].first << endl;
		}
		else
		{
			char c[MAXRSIZE];
			size_t len = fread(&c, 1, MAXRSIZE, fd), i = 0;
			uint64_t _km_f = 0, _km_r = 0;
			bool _isfull = false;
			uint8_t cpt = 0;
			if (len > 0 && c[0] == '>')
			{
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							if (_isfull)
							{
								_km_f >>= 2;
								_km_f += m_powerTable[cpt][m_table[c[i]]];

								/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_f, m_targetsID[t].second, 1);
								///////////////////////////////////////////////////////////////////////
								i++;
								continue;
							}
							_km_r <<= 2;    
							_km_r += 3 - m_table[c[i]];
							if (cpt == m_kmerSize - 1)
							{       _isfull = true;
								//[/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_r, m_targetsID[t].second, 1);
								//]//////////////////////////////////////////////////////////////////////
								_km_f = _km_r;
								// The following 6 lines come from Jellyfish source code
								_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
								_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
								_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
								_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
								_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
								_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
								i++;
								continue;
							}
							cpt++;  i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{       i++;
							continue;   }
						if (m_table[c[i]] == -1)
						{        nt++; _km_r = 0; cpt = 0; _isfull = false; i++;
							continue;  }
						if (m_table[c[i]] == -2)
						{
							_km_r = 0; cpt = 0; _isfull = false;
							while (len != 0)
							{
								while(i < len && c[i] != '\n')  {       i++;}
								if (i < len)            {   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")    ";
				continue;
			}
			if (len > 0 && c[0] == '@')
			{
				while (len != 0)
				{
					while(i < len && c[i] != '\n')            {       i++;}
					if (i < len)      {   break;}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				i++;
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							if (_isfull)
							{
								_km_f >>= 2; 
								_km_f += m_powerTable[cpt][m_table[c[i]]];

								/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_f, m_targetsID[t].second, 1);
								///////////////////////////////////////////////////////////////////////
								i++;
								continue;
							}
							_km_r <<= 2;    
							_km_r += 3 - m_table[c[i]];
							if (cpt  == m_kmerSize - 1)
							{       _isfull = true;
								//[/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_r, m_targetsID[t].second, 1);
								//]//////////////////////////////////////////////////////////////////////
								_km_f = _km_r;
								// The following 6 lines come from Jellyfish source code
								_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
								_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
								_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
								_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
								_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
								_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));

								i++;
								continue;
							}
							cpt++;
							i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{
							_km_r = 0; cpt = 0; _isfull = false;
							// Pass third line
							i++;
							while (len != 0)
							{
								while(i < len && c[i] != '\n')          {       i++;}
								if (i < len)            {   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass fouth line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')         {       i++;}
								if (i < len)         {   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass first line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')        {       i++;}
								if (i < len){   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						if (m_table[c[i]] == -1)
						{       nt++;  _km_r = 0; cpt = 0; _isfull = false;
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
				continue;
			}
			// spectrum form
			fseek(fd, 0, 0);
			string s_kmer = "";
			ITYPE val;
			while (getFirstAndSecondElementInLine(fd, s_kmer, val))
			{       if (val > m_minCountTarget)
				{       commonKmersHT.addElement(s_kmer, m_targetsID[t].second, (size_t) val);}
			}
			fclose(fd);
			cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
		}
	}
	cerr << nt << " nt read in total." << endl;
	cerr << "Mother Hashtable successfully built. " << commonKmersHT.Size() << " " << m_kmerSize << "-mers stored." <<  endl;
	size_t sizeMotherTable = commonKmersHT.Size();
	//
	commonKmersHT.SaveIntersectionMultiple(_filesHTC, m_labels_c);
	commonKmersHT.SaveMultiple(_filesHT, m_labels);
	//
	commonKmersHT.SortAllHashTable(2);
	commonKmersHT.RemoveCommon(m_labels_c, m_minCountTarget);
	char * cfname = (char*) calloc(130, sizeof(char));
	getdbName(cfname);
	cerr << "Creating database in disk..." << endl;
	uint64_t nbElement = commonKmersHT.Write(cfname,2);
	free(cfname);
	cfname = NULL;
	cerr << nbElement << " " << m_kmerSize << "-mers successfully stored in database." << endl;

	return sizeMotherTable;
	///////////////////////////////////////////////////////////////////////////////
}

/**
 * Process input file and classify its sequences.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::getObjectsDataComputeFullGPU(const uint8_t * _map,  const size_t&   nb, const char* _fileResult)
{
	// find and store the start and end of each reads name,
	// the start and end of each read and its length
	size_t i_r = 0, bigSteps = 1 + nb / m_numBatches, lastSize = 0;
	if (_map[0] == '>')			// fasta
	{
#ifdef _OPENMP
#pragma omp parallel for private(i_r) firstprivate(lastSize)
#endif
		for (i_r = 0; i_r < m_numBatches ; i_r++)
		{
			size_t i = bigSteps * i_r, i_l = 0, bigMax = bigSteps*(i_r+1);
			
			vector<size_t>	readsLength;
			vector<size_t>	seqSNames;
			vector<size_t>	seqENames;
			vector<size_t>	readsSPos;
			vector<size_t>	readsEPos;
			
			if (lastSize)
			{
				readsLength.reserve(lastSize);
				seqSNames.reserve(lastSize);
				seqENames.reserve(lastSize);
				readsSPos.reserve(lastSize);
				readsEPos.reserve(lastSize);
			}	
				
			while (_map[i++] != '>')
			{}
			while (true)
			{
				readsLength.push_back(0);
				seqSNames.push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				seqENames.push_back(i);

				while (i < nb && _map[i++] != '\n')
				{}

				readsSPos.push_back(i);
				readsEPos.push_back(i);
				while (i < nb && _map[i]!= '>')
				{
					while (i < nb && _map[i]!= '\n')
					{
						i++;
					}
					readsLength[i_l]++;
					readsEPos[i_l] = i++;
				}
				readsLength[i_l] *= -1;
				readsLength[i_l] += readsEPos[i_l] -  readsSPos[i_l] + 1;
				if ((i >= bigMax && _map[i] == '>' )|| i >= nb)
				{       break;}
				i++;
				i_l++;
			}
			
			m_readsLength[i_r].swap(readsLength);
			m_seqSNames[i_r].swap(seqSNames);
			m_seqENames[i_r].swap(seqENames);
			m_readsSPos[i_r].swap(readsSPos);
			m_readsEPos[i_r].swap(readsEPos);
			
			lastSize = m_readsLength[i_r].size()*1.05;
		}
	}
	else if (_map[0] == '@')		// fastq
	{
		m_posReads[0] = 1;
		size_t pos1, pos2, pos3, pos4, pos5, pos6;
		for(i_r = 1; i_r < m_numBatches ; i_r++)
		{
			size_t i = bigSteps * i_r;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos1 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos2 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos3 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos4 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos5 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos6 = i;
			if (_map[pos1] == '@')
			{
				i = pos2;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos3 && _map[i] == '+')
				{
					m_posReads[i_r] = pos1+1;
					continue;
				}
			}
			if (_map[pos2] == '@')
			{
				i = pos3;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos4 && _map[i] == '+')
				{       m_posReads[i_r] = pos2+1;
					continue;
				}
			}
			if (_map[pos3] == '@')
			{
				i = pos4;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos5 && _map[i] == '+')
				{       m_posReads[i_r] = pos3+1;
					continue;
				}
			}
			if (_map[pos4] == '@')
			{
				i = pos5;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos6 && _map[i] == '+')
				{       m_posReads[i_r] = pos4+1;
					continue;
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for private(i_r) firstprivate(lastSize)
#endif
		for (i_r = 0; i_r < m_numBatches ; i_r++)
		{
			size_t iNext = i_r+1 < m_numBatches ? m_posReads[i_r+1]: nb;
			size_t i = m_posReads[i_r], i_l = 0;

			vector<size_t>	readsLength;
			vector<size_t>	seqSNames;
			vector<size_t>	seqENames;
			vector<size_t>	readsSPos;
			vector<size_t>	readsEPos;
				
			if (lastSize)
			{
				readsLength.reserve(lastSize);
				seqSNames.reserve(lastSize);
				seqENames.reserve(lastSize);
				readsSPos.reserve(lastSize);
				readsEPos.reserve(lastSize);
			}

			while (true)
			{
				seqSNames.push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				seqENames.push_back(i);

				while (i < nb && _map[i++] != '\n')
				{}

				readsSPos.push_back(i);
				readsEPos.push_back(i);
				while (i < nb && _map[i]!= '\n')
				{
					i++;
				}
				readsEPos[i_l] = i++;
				readsLength.push_back(readsEPos[i_l] -  readsSPos[i_l]);
				// Pass third line
				while (i < nb && _map[i++] != '\n')
				{}
				// Pass fourth line
				while (i < nb && _map[i++] != '\n')
				{}

				if ((++i) >= iNext)
				{       break;}
				i_l++;
			}
			
			m_readsLength[i_r].swap(readsLength);
			m_seqSNames[i_r].swap(seqSNames);
			m_seqENames[i_r].swap(seqENames);
			m_readsSPos[i_r].swap(readsSPos);
			m_readsEPos[i_r].swap(readsEPos);
			
			lastSize = m_readsLength[i_r].size()*1.05;
		}
	}
	else 
	{ 
		cerr << "Failed to recognize the format of the file." << endl; exit(-1) ;
	}
	
	// get index of first read in batch, count the number of all reads
	vector<ITYPE> indexBatches;
	indexBatches.resize(m_numBatches+1);
	indexBatches[0] = 0;
	size_t maxReads = 0;
	for(size_t i = 0; i < m_numBatches; i++)
	{
		//~ indexReads[i].reserve(m_readsLength[i].size());
		indexBatches[i+1] = indexBatches[i] + m_readsLength[i].size();
		if(m_readsLength[i].size() > maxReads) maxReads = m_readsLength[i].size();
	}
	m_nbObjects = indexBatches.back();
	
	// get average read length for memory allocation
	size_t numContainerMax = 0;
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for(size_t i = 0; i < m_numBatches; i++)
	{
		size_t avgReadLength = 0;
		for(size_t n = 0; n < m_readsLength[i].size(); n++)
		{
			avgReadLength += m_readsLength[i][n];
		}	
		avgReadLength = avgReadLength/m_readsLength[i].size() +1;
		size_t numContainer = avgReadLength/(sizeof(CONTAINER)*4)+3;
#ifdef DEBUG_BATCH
		cerr << "Batch " << i << ": AVG read length " << avgReadLength
			<< ", Estimated # containers per read: " << numContainer << "\n";
#endif			

		if(numContainer > numContainerMax)
		{
#ifdef _OPENMP
			#pragma omp critical (numContainer)
			if(numContainer > numContainerMax)
#endif			
				numContainerMax = numContainer;
		}
	}
#ifdef DEBUG_BATCH
	cerr << "Estimated # containers per read: " << numContainerMax << "\n";
#endif			
	
	// compact results
	// need to store index and counter for each target hit, also sum
	if (m_isLightLoading)
		m_resultRowSize = 2*MAXHITS+2; 		// = 48, maximum used in light test 10M was 1+21*2=43
	else
		m_resultRowSize = 2*MAXHITS+2; 		// = 32, maximum used in normal test 10M was 1+2*10 = 21
		
	// space for result vector (sumN, indexBest, best, index_sBest, s_best)
	m_finalResultsRowSize = 5;
	
	// the pointers to the containers
	std::vector<uint32_t*> readsPointer;
	// the part lengths and nucleotides for each read
	std::vector<CONTAINER*> readsInContainers;
	// allocate gpu mem
	size_t total = m_cuClarkDb->malloc(m_nbObjects, maxReads, maxReads*numContainerMax,
										indexBatches,
										m_fullResults, m_resultRowSize,
										m_finalResults, m_finalResultsRowSize,
										m_isExtended,
										readsPointer,
										readsInContainers);
	
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		
#ifdef _OPENMP
		#pragma omp for schedule(dynamic) //ordered//schedule(static,1) //private(i_r)
#endif
		for (i_r = 0; i_r < m_numBatches; i_r++)
		{
			CONTAINER _kmerContainer;
			// going to store 4 nucleotides per byte
			size_t nucsPerContainer = sizeof(_kmerContainer) *4;
			size_t partBegin;
			size_t curNucs = 0;
			bool _newSeq = true;
			size_t i_lr = 0, i_c = 0, i_s = 0;	
			
			uint32_t* myReadsPointer = readsPointer[i_r];
			CONTAINER* myReadsInContainers = readsInContainers[i_r];
			size_t containerCount = 0;
			
			while (i_lr < m_readsLength[i_r].size())
			{
				// check if read has no kmers
				i_c = m_readsLength[i_r][i_lr] < m_kmerSize ? m_readsEPos[i_r][i_lr]: m_readsSPos[i_r][i_lr];
				
				partBegin = containerCount;
				myReadsPointer[i_lr] = partBegin;
				// store read in containers
				while (i_c < m_readsEPos[i_r][i_lr])
				{
					// part continues
					if (m_table[_map[i_c]] >= 0)
					{
						if (_newSeq)
						{
							// add element to store part length
							myReadsInContainers[containerCount++] = 0;		
							if (myReadsInContainers[partBegin] < m_kmerSize)
							{
								containerCount = partBegin+1;
								myReadsInContainers[partBegin] = 0;
							}
							else
							{
								// advance to new seq
								partBegin = containerCount-1;
							}
							
							
							_newSeq = false;
						}
						_kmerContainer <<= 2;
						_kmerContainer ^= m_rTable[_map[i_c]];
						curNucs++;
						i_c++;
						// store container if full and update part length
						if ( curNucs == nucsPerContainer )
						{
							myReadsInContainers[containerCount++] = _kmerContainer;
							myReadsInContainers[partBegin] += curNucs;
							curNucs = 0;
						}
						continue;
					}
					if (_map[i_c] == '\n')
					{
						i_c++;
						continue;
					}
					// part ends
					// store container and update part length
					if (curNucs > 0)
					{
						_kmerContainer <<= 2*(nucsPerContainer-curNucs);
						myReadsInContainers[containerCount++] = _kmerContainer;
						myReadsInContainers[partBegin] += curNucs;
					}
					// reset for new seq
					_kmerContainer = 0; curNucs = 0; _newSeq = true;
					i_c++;
				}
				// read & part end
				// store container and update part length
				if (curNucs > 0)
				{
					_kmerContainer <<= 2*(nucsPerContainer-curNucs);
					myReadsInContainers[containerCount++] = _kmerContainer;
					myReadsInContainers[partBegin] += curNucs;
				}
				// remove seq if it's too short
				// this happens if read length >= kmer size but all seq in read < kmer size
				if (myReadsInContainers[partBegin] < m_kmerSize)
				{
					containerCount = partBegin;
				}
				// reset for new read
				i_lr++;
				_kmerContainer = 0; curNucs = 0; _newSeq = true;	
			}
			
			if (containerCount > (uint32_t)-1)
			{
				cerr << "ERROR: Batch overflow. Please increase the number of batches (-b <numberofbatches>)." << endl;
				exit(-1);
			}	
			
			myReadsPointer[m_readsLength[i_r].size()] = containerCount;
			
			float numContainerActual = (float)containerCount/m_readsLength[i_r].size();			
#ifdef DEBUG_BATCH
			cerr << "Batch " << i_r << "\t # Containers: " << containerCount
				 << " (AVG # per read: " << numContainerActual << ").\n";
#endif			 
			if (numContainerActual > numContainerMax)
			{
				cerr << "Bad container estimation. Abort." << endl;
				exit(-1);
			}		
#ifdef DEBUG_BATCH
			cerr 	<< "Host " << omp_get_thread_num()
					<< " Batch " << i_r << "\t CUDA start.\t"
					<< "Read data size: " << readsInContainers.size()*sizeof(CONTAINER) /1000 /1000.0 << " MB\t"
					<< "Result data size: " << m_finalResultsRowSize*m_readsLength[i_r].size() /1000 /1000.0 << " MB\n";
#endif				
			// pass batch parameters to m_cuClarkDb			
			m_cuClarkDb->readyBatch(i_r, m_readsLength[i_r].size(), containerCount);

			// query batches sequencially
#ifdef _OPENMP
			#pragma omp critical(cudaQuery)
			//~ #pragma omp ordered
#endif
			{
				m_batchScheduled[i_r] = m_cuClarkDb->queryBatch( i_r, m_isExtended);
				
#ifdef DEBUG_BATCH
				// print batch information
				size_t batchSize = (m_readsLength[i_r].size()+1)*sizeof(uint32_t) + containerCount*sizeof(CONTAINER);
				cerr << "Batch " << i_r << " scheduled." 
						  << " Objects: " << m_readsLength[i_r].size()
						  << " Size: " << batchSize/1000/1000.0 << " MB"
						  << endl;
#endif
			}
					
			// start printing results early if possible
#ifdef _OPENMP
			if (m_dbParts == 1 && omp_get_max_threads() > 1 && i_r == 0)
			{
				printExtendedResultsSynced(_map, _fileResult);
			}
#endif
			
		}
	}

	// swap db parts if available (multi threaded)
	while(m_cuClarkDb->swapDbParts())
	{			
		// query batches again
		for (i_r = 0; i_r < m_numBatches; i_r++)
		{
			m_batchScheduled[i_r] = m_cuClarkDb->queryBatch( i_r, m_isExtended, true);
		}
	}
	
	// print results
#ifdef _OPENMP
	if (m_dbParts > 1 || omp_get_max_threads() == 1)
#endif
	{
		//~ std::cerr << "CPU ready to write.\n";
		printExtendedResultsSynced(_map, _fileResult);
	}
	
	return;
}

/**
 * Get targets data from targets files.
 * Return true if database is already available.
 */
template <typename HKMERr>
bool CuCLARK<HKMERr>::getTargetsData(const char* _filesName, vector<string>& _filesHT, vector<string>& _filesHTC, const bool& _creatingkmfiles, const ITYPE& _samplingfactor)
{
	FILE * meta_f = fopen(_filesName,"r");
	if (meta_f == NULL)
	{
		cerr << "Failed to open targets data in file: " << _filesName << endl;
		exit(-1);
	}

	string subfile = "";
	while (getFirstElementInLineFromFile(meta_f, subfile))
	{
		FILE * t_f = fopen(subfile.c_str(),"r");
		if (t_f == NULL)
		{
			cerr << "Failed to open file: " << subfile << " defined in " << _filesName << endl;
			exit(-1);
		}
		fclose(t_f);
	}
	fclose(meta_f);	
	meta_f = fopen(_filesName,"r");

	bool areHTfilespresent = true;
	while(getLineFromFile(meta_f, subfile))
	{
		vector<string> ele;
		getElementsFromLine(subfile, 3, ele);
		pair< string, string > target;

		target.first = ele[0];

		if (ele.size() > 1)
		{
			target.second = ele[1];
			std::vector<string>::iterator it = find (m_labels.begin(), m_labels.end(), ele[1]); 
			if ( it  == m_labels.end())
			{	
				m_labels.push_back(ele[1]); 
			}
			char * fname = (char*) calloc(100, sizeof(char));
			sprintf(fname, "%s/%s_k%lu.ht", m_folder, ele[1].c_str(), (size_t)m_kmerSize);
			FILE * fd = fopen(fname, "r");
			areHTfilespresent =  (fd != NULL) && areHTfilespresent;
			free(fname);
			fname = NULL;
			if (fd != NULL)
			{
				fclose(fd);
			}
		}
		else
		{	cerr << " Missing label for " << ele[0] <<  endl; exit(-1); }

		if (ele.size() > 2)
		{
			std::vector<string>::iterator it = find (m_labels_c.begin(), m_labels_c.end(), ele[2]);
			if ( it  == m_labels_c.end())
			{       m_labels_c.push_back(ele[2]);} 

			char * fname = (char*) calloc(100, sizeof(char));
			sprintf(fname, "%s/%s_k%lu.ht", m_folder, ele[2].c_str(), (size_t) m_kmerSize);
			FILE * fd = fopen(fname, "r");
			areHTfilespresent =  (fd != NULL) && areHTfilespresent;
			free(fname);
			fname = NULL;
			if (fd != NULL)
			{
				fclose(fd);
			}

		}
		m_targetsID.push_back(target);
	}
	fclose(meta_f);
	print(_creatingkmfiles, _samplingfactor);

	if (_creatingkmfiles)
	{
		createTargetFilesNames(_filesHT, _filesHTC);
	}

	size_t sHTS = m_labels_c.size() + m_labels.size();

	m_targetsName.push_back("NA");
	for(size_t t=0 ; t < m_labels.size(); t++)
	{       m_targetsName.push_back(m_labels[t]);       }
	for(size_t t=0 ; t < m_labels_c.size(); t++)
	{       m_targetsName.push_back(m_labels_c[t]);       }

	char * cfname = (char*) calloc(130, sizeof(char));
	char * cfname_s = (char*) calloc(130+4, sizeof(char));
	getdbName(cfname);
	sprintf(cfname_s, "%s.sz", cfname);
	FILE * dbfd_sz = fopen(cfname_s, "r");
	sprintf(cfname_s, "%s.ky", cfname);
	FILE * dbfd_ky = fopen(cfname_s, "r");
	sprintf(cfname_s, "%s.lb", cfname);
	FILE * dbfd_lbl = fopen(cfname_s, "r");
	if (dbfd_sz != NULL && dbfd_ky != NULL && dbfd_lbl != NULL)
	{
		areHTfilespresent = true;
		fclose(dbfd_sz);
		fclose(dbfd_ky);
		fclose(dbfd_lbl);
	}
	free(cfname);
	free(cfname_s);
	cfname = NULL;
	cfname_s = NULL;
	return areHTfilespresent;
}

/**
 * Prints program intro to standard output, additional parameters to standard error.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::print(const bool& _creatingkmfiles, const ITYPE& _samplingfactor) const
{
	cerr << "CuCLARK version " << VERSION << " (Copyright 2016 Robin Kobus, rkobus@students.uni-mainz.de)" << endl;
	cerr << "Based on CLARK version 1.1.3 (UCR CS&E. Copyright 2013-2016 Rachid Ounit, rouni001@cs.ucr.edu) " << endl;
	if (m_minCountTarget > 0)
	{
		cerr << "Minimum k-mers occurences in Targets is set to " << m_minCountTarget << endl;
	}
	if (_creatingkmfiles)
	{
		cerr << "Creation of targets specific k-mers files requested " << endl;
	}
	if (m_isLightLoading)
	{
		cerr << "Using light database in RAM (" << m_iterKmers << ")" << endl;
	}
	if (_samplingfactor > 2)
	{
		cerr << "Sampling factor is " << _samplingfactor << endl;
	}
}

/**
 * Prints speed stats for processing an input file to standard output.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::printSpeedStats(const struct timeval& _requestEnd, const struct timeval& _requestStart, const char* _fileResult) const 
{
	double diff = (_requestEnd.tv_sec - _requestStart.tv_sec) + (_requestEnd.tv_usec - _requestStart.tv_usec) / 1000000.0;
	cout <<" - Assignment time: "<<diff<<" s. Speed: ";
	cout << (size_t) (((double) m_nbObjects)/(diff)*60.0)<<" objects/min. ("<< m_nbObjects<<" objects)."<<endl;
	cout <<" - Results stored in " << _fileResult << endl;
}

/**
 * Prints results to file in normal or extended format.
 * Waits for a batch to finisch before printing its results.
 */
template <typename HKMERr>
void CuCLARK<HKMERr>::printExtendedResultsSynced(const uint8_t * _map,  const char* _fileResult)
{
	ofstream f_out;
	f_out.open(_fileResult, std::ofstream::out );

	// print header
	string header[] = {"Length", "Gamma","1st_assignment", "score1", "2nd_assignment", "score2", "confidence"};
	size_t headerSize = 7;

	f_out <<"Object_ID";

	if (m_isExtended)
	{
		for(size_t t = 1 ; t < m_targetsName.size();  t++)
		{
			f_out << "," << m_targetsName[t];
		}
	}

	for(size_t t = 0 ; t < headerSize ;  t++)
	{       f_out << "," << header[t]; }
	f_out << endl;

	f_out.close();
	//
	FILE *fout = fopen(_fileResult, "a");

	ITYPE best = 0, s_best = 0, indexBest = 0, index_sBest = 0, total = 0;
	double gamma = 0, delta = 0;
	
	//// print extended results
	if (m_isExtended)
	{
		int nonzero_count = 0, nonzero_sum = 0,
			nonzero_max = 0, nonzero_min = m_targetsName.size()-1;
		int targetIndex, targetScore;
		size_t writeIndex;
		
		int i_r=0;
		int i_lr=0;
		char objectName[OBJECTNAMEMAX];
		ITYPE objectNorm;
		size_t batchLastIndex = m_readsLength[i_r].size();
		
		// wait for first batch
		while(!m_batchScheduled[i_r]);
		m_cuClarkDb->waitForBatch(i_r);
		cerr << "Writing extended results... " << endl;
		
		for(size_t t = 0; t < m_nbObjects; t++)
		{
			if(t==batchLastIndex)
			{
				// wait for next batch
				while(!m_batchScheduled[++i_r]);
				m_cuClarkDb->waitForBatch(i_r);
				batchLastIndex += m_readsLength[i_r].size();
				i_lr=0;
			}
			
			writeIndex = 0;
			std::ostringstream ss;
			
			// print all scores
			for(size_t r_h = 0 ; r_h < m_fullResults[t*m_resultRowSize] ; r_h++)
			{
				// get next target
				targetIndex = m_fullResults[t*m_resultRowSize+r_h*2+1];
				targetScore = m_fullResults[t*m_resultRowSize+r_h*2+2];

				// write zeros until target
				for ( ; writeIndex < targetIndex; writeIndex++)
					ss << ",0";
				// write target score
				ss << "," << targetScore;
				writeIndex++;
					
			}
			// write zeros until end
			for ( ; writeIndex < m_targetsName.size()-1; writeIndex++)
				ss << ",0";

			total 		= m_finalResults[t*m_finalResultsRowSize];
			indexBest 	= m_finalResults[t*m_finalResultsRowSize+1];
			best 		= m_finalResults[t*m_finalResultsRowSize+2];
			index_sBest = m_finalResults[t*m_finalResultsRowSize+3];
			s_best 		= m_finalResults[t*m_finalResultsRowSize+4];

			size_t nameSize = m_seqENames[i_r][i_lr] - m_seqSNames[i_r][i_lr];
			if (nameSize >= OBJECTNAMEMAX) nameSize = OBJECTNAMEMAX-1;
			strncpy(objectName, (const char*) &_map[m_seqSNames[i_r][i_lr]], nameSize);
			objectName[nameSize] = '\0';
		
			objectNorm = m_isPaired ? m_readsLength[i_r][i_lr] - NBN : m_readsLength[i_r][i_lr];
			i_lr++;
			
// remove 2. double cast (not done to keep comparability)
// also float would be totally sufficient...
// same applies below
			//~ float gamma  = (float)(total)/((objectNorm - m_kmerSize) + 1.0);
			//~ float delta = best + s_best;
			//~ delta = (delta < 0.001) ? 0: ((float) best)/(delta);

			gamma = (double)total / ((double)objectNorm - m_kmerSize + 1.0);
			delta = best + s_best;
			delta = (delta < 0.001) ? 0: (double) best/ delta;
			
			
			// print name, hits, length, hit rate, best, second best, confidence score
			fprintf(fout,"%s%s,%u,%g,%s,%u,%s,%u,%g\n",
					objectName,ss.str().c_str(),objectNorm,gamma,
					m_targetsName[indexBest].c_str(),best,m_targetsName[index_sBest].c_str(),s_best,
					delta);
					
			/// extra target info		
			nonzero_count = m_fullResults[t*m_resultRowSize];
			if (nonzero_count > nonzero_max) nonzero_max = nonzero_count;
			if (nonzero_count < nonzero_min) nonzero_min = nonzero_count;
			nonzero_sum += nonzero_count;
			///
		}
		fclose(fout);
		cerr << "Done." << endl;
		
		/// extra target info
		cerr << "MIN targets: " << nonzero_min
			 << ", MAX targets: " << nonzero_max
			 << ", AVG targets: " << (float)nonzero_sum / m_nbObjects
			 << "\n";
		///

		return;
	}
	
	//// print normal results
	int i_r=0;
	int i_lr=0;
	char objectName[OBJECTNAMEMAX];
	ITYPE objectNorm;
	size_t batchLastIndex = m_readsLength[i_r].size();
	
	// wait for first batch
	while(!m_batchScheduled[i_r]);
	m_cuClarkDb->waitForBatch(i_r);
	cerr << "Writing results... " << endl;
	
	for(size_t t = 0; t < m_nbObjects; t++)
	{
		if(t==batchLastIndex)
		{
			// wait for next batch
			while(!m_batchScheduled[++i_r]);
			m_cuClarkDb->waitForBatch(i_r);
			batchLastIndex += m_readsLength[i_r].size();
			i_lr=0;
		}
			
		total 		= m_finalResults[t*m_finalResultsRowSize];
		indexBest 	= m_finalResults[t*m_finalResultsRowSize+1];
		best 		= m_finalResults[t*m_finalResultsRowSize+2];
		index_sBest = m_finalResults[t*m_finalResultsRowSize+3];
		s_best 		= m_finalResults[t*m_finalResultsRowSize+4];
		
		size_t nameSize = m_seqENames[i_r][i_lr] - m_seqSNames[i_r][i_lr];
		if (nameSize >= OBJECTNAMEMAX) nameSize = OBJECTNAMEMAX-1;
		strncpy(objectName, (const char*) &_map[m_seqSNames[i_r][i_lr]], nameSize);
		objectName[nameSize] = '\0';
		
		objectNorm = m_isPaired ? m_readsLength[i_r][i_lr] - NBN : m_readsLength[i_r][i_lr];
		i_lr++;
		
// see comments above
		//~ float gamma  = (float)(total)/((objectNorm - m_kmerSize) + 1.0);
		//~ float delta = best + s_best;
		//~ delta = (delta < 0.001) ? 0: ((float) best)/(delta);
		
		gamma = (double)(total)/(((double) objectNorm - m_kmerSize) + 1.0);
		delta = best + s_best;
		delta = (delta < 0.001) ? 0: ((double) best)/(delta);			
		
		// print name, length, hit rate, best, second best, confidence score
		fprintf(fout,"%s,%u,%g,%s,%u,%s,%u,%g\n",
				objectName,objectNorm,gamma,
				m_targetsName[indexBest].c_str(),best,m_targetsName[index_sBest].c_str(),s_best,
				delta);
	}
	fclose(fout);
	cerr << "Done." << endl;	
}
