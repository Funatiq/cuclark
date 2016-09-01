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
 * New file:
 * New database structure and GPU handling for database and queries.
 * Includes kernels to query batches of reads, merge results and find best results.
 */
 
#include "CuClarkDB.cuh"

#include <iostream>
#include <fstream>
#include <iterator>
//~ #include <numeric>	// partial_sum
#include <cstring>	// memcpy

// for debugging prints
//~ #include <inttypes.h>	// PRIu64 print makro
//~ #include <bitset>		// print kmer container

#define CUERR {														\
	cudaError_t err;												\
	if ((err = cudaGetLastError()) != cudaSuccess)					\
	{																\
		std::cerr << "CUERR '" << cudaGetErrorString(err) << "' in "\
				  << __FILE__ << ", line " << __LINE__ << "\n";		\
		exit(1);													\
	}																\
}

#define CUMEMERR {													\
	if (cudaGetLastError()== cudaErrorMemoryAllocation)				\
	{																\
		std::cerr << "ERROR: Out of GPU memory.\n"					\
				  << "Please increase the number of batches "		\
				  << "(-b <numberofbatches>).\n";					\
		exit(1);													\
	}																\
}

// forward declaration
template <typename HKMERr>
__device__ bool queryElement (const uint8_t& k, const uint64_t& _ikmer,
			uint32_t* d_bucketPointers, HKMERr* d_keys, ILBL* d_labels,
			uint32_t dbPartStart, uint32_t dbPartEnd,
			//~ uint8_t dbParts, uint8_t dbPart,
			ILBL& _returnLabel);
template <typename HKMERr>
__global__ void queryKernel (uint8_t k,
			uint32_t* readsPointer, CONTAINER* readsInContainers,
			uint32_t* bucketPointers, HKMERr* keys, ILBL* labels,
			uint32_t dbPartStart, uint32_t dbPartEnd,
			//~ uint8_t dbParts, uint8_t dbPart,
			RESULTS* results, size_t pitch, size_t numTargets);
__global__ void mergeKernel (RESULTS* resultA, RESULTS* resultB, size_t pitch, size_t numReads, RESULTS* results);
__global__ void resultKernel (RESULTS* scores, size_t spitch, size_t numReads, RESULTS* results, size_t rpitch);

/**
 * Constructor:
 * Initialize variables, find CUDA devices
 */	
template <typename HKMERr>
CuClarkDB<HKMERr>::CuClarkDB(const size_t _numDevices, const uint8_t _k, const size_t _numBatches, const size_t _numTargets)
							: m_k(_k),m_numTargets(_numTargets),m_numBatches(_numBatches)
{	
	m_numReads.resize(m_numBatches);
	m_sizeReadsPointer.resize(m_numBatches);
	m_sizeReadsInContainers.resize(m_numBatches);
	
	h_readsPointer.resize(m_numBatches);
	h_readsInContainers.resize(m_numBatches);
	
	h_results.resize(m_numBatches);
	h_resultsFinal.resize(m_numBatches);
	
	m_batchFinishedEvents.resize(m_numBatches);
	for(int i=0; i <m_numBatches; i++)
	{
		cudaEventCreateWithFlags(&m_batchFinishedEvents[i],cudaEventDisableTiming);
	}
	
	m_numDevices = 0;
	m_dbPartsPerDevice = 1;
	
	// cf. CUDA samples/0_Simple/simpleP2P
	std::cerr << "Checking for CUDA devices: ";
	cudaGetDeviceCount(&m_numDevices);
	CUERR
	if (m_numDevices > 0)
		std::cerr << m_numDevices << " device(s) found.\n";
	else
	{
		std::cerr << "No CUDA devices found. Abort.\n";
		exit(1);
	}
	
	std::vector<cudaDeviceProp> prop(m_numDevices);
	for(int i=0; i<m_numDevices; i++)
	{
		cudaGetDeviceProperties(&prop[i], i);
		CUERR
		std::cerr << "Device " << i << " = " << prop[i].name << "\n";
	}
	
	if (m_numDevices < _numDevices)
	{
		std::cerr << _numDevices << " CUDA devices requested. Insufficient devices found. Abort.\n";
		exit(1);
	}
	
	if (_numDevices > 0)
	{
		std::cerr << "Using " << _numDevices << " CUDA devices as requested.\n";
		m_numDevices = _numDevices;
	}
	
	for(int i=0; i<m_numDevices; i++)
	{
		if(strcmp(prop[i].name,"GeForce GTX TITAN X") == 0)
			m_dbPartsPerDevice = DBPARTSPERDEVICE;
	}
	
	std::vector<int>	gpuid(m_numDevices);
	int gpu_count = 0;
	m_memSizes.resize(m_numDevices*m_dbPartsPerDevice);
	
	for(int i=0; i<m_numDevices; i++)
	{
		// locate devices capable of Peer-to-Peer
		if(prop[i].major >= 2)
		{
			gpuid[gpu_count++] = i;
			//~ std::cerr << "Device " << i << " = " << prop[i].name << " is capable of P2P.\n";
			//~ printf("> GPU%d = \"%15s\" is capable of Peer-to-Peer (P2P)\n", i, prop[i].name);
		}
		

		size_t freeMem, totalMem;
		cudaSetDevice(i);
		cudaMemGetInfo(&freeMem, &totalMem);
		CUERR
		//~ std::cerr << "Device " << i << " free: " << freeMem/1000000
									//~ << " total: " << totalMem/1000000
									//~ << " global: " << prop[i].totalGlobalMem/1000000 
									//~ << std::endl;
		if (freeMem < 1000000000)
		{
			std::cerr << "Device " << i << " has less than 1GB of free memory. Abort.\n";
			exit(1);
		}
		
		freeMem -= RESERVED;
		for (int j=0; j<m_dbPartsPerDevice; ++j)
		{
			m_memSizes[m_dbPartsPerDevice*i+j] = freeMem/m_dbPartsPerDevice + (i+j>0 ? m_memSizes[m_dbPartsPerDevice*i+j-1] : 0);		
		}
	}
	
	int can_access_peer;
	// check all the combinations of supported P2P GPUs
    for (int i = 0; i < gpu_count; i++)
    {
        for (int j = i+1; j < gpu_count; j++)
        {
            if (gpuid[i] == gpuid[j])
            {
                continue;
            }
            cudaDeviceCanAccessPeer(&can_access_peer, gpuid[i], gpuid[j]);
            CUERR
            
            if(can_access_peer)
            {
				std::cerr << "Enabling peer access between devices " << gpuid[i] << " and " << gpuid[j] << "\n";
				cudaSetDevice(gpuid[i]);
				cudaDeviceEnablePeerAccess(gpuid[j], 0);
				CUERR
				cudaSetDevice(gpuid[j]);
				cudaDeviceEnablePeerAccess(gpuid[i], 0);
				CUERR
			}
        }
    }
    
    // pointers for each device
    d_bucketPointers.resize(m_numDevices*m_dbPartsPerDevice);
	d_keys.resize(m_numDevices*m_dbPartsPerDevice);
	d_labels.resize(m_numDevices*m_dbPartsPerDevice);
	
	d_readsPointer.resize(m_numDevices);
	d_readsInContainers.resize(m_numDevices);
	
	d_results.resize(m_numDevices);

	for (int i=0; i<m_numDevices; i++)
	{
		if (m_dbPartsPerDevice > 1)
		{// results for each part, +1 for merging
			d_results[i].resize(m_dbPartsPerDevice+1);
		}
		else
		{// 3 results on every other device for merging
			if (i%2 == 0)
				d_results[i].resize(3);
			else
				d_results[i].resize(1);
		}		
	}
	
	d_resultsFinal.resize(m_numDevices);
}

/**
 * Destructor:
 * Free memory allocations
 */
template <typename HKMERr>
CuClarkDB<HKMERr>::~CuClarkDB() 
{
	for (int i=0; i<m_dbParts; i++)
	{		
		cudaFreeHost(h_bucketPointers[i]);
		cudaFreeHost(h_keys[i]);
		cudaFreeHost(h_labels[i]);
	}

	for(int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		
		for(int j=0; j<m_dbPartsPerDevice; ++j)
		{
			int index = m_dbPartsPerDevice*i+j;
			cudaFree(d_bucketPointers[index]);
			cudaFree(d_keys[index]);
			cudaFree(d_labels[index]);
			//~ CUERR
		}
	}
	
	for(int i=0; i<m_batchFinishedEvents.size(); i++)
		cudaEventDestroy(m_batchFinishedEvents[i]);
	//~ CUERR
	
	for(int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		CUERR
		cudaDeviceReset();
	}
}

/**
 * Reset memory allocations before new input file.
 */
template <typename HKMERr>
void CuClarkDB<HKMERr>::free()
{
	//~ cudaSetDevice(0);
	for (int i=0; i<m_numBatches; i++)
	{
		cudaFreeHost(h_readsPointer[i]);
		cudaFreeHost(h_readsInContainers[i]);
		CUERR
	}
	
	cudaFreeHost(h_results[0]);
	cudaFreeHost(h_resultsFinal[0]);
	CUERR	
	
	for(int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
		CUERR
		
		cudaFree(d_readsPointer[i]);
		cudaFree(d_readsInContainers[i]);
		CUERR
		
		for (int j=0; j<d_results[i].size(); j++)
			cudaFree(d_results[i][j]);
		CUERR
		
		cudaFree(d_resultsFinal[i]);
		CUERR
	}

	//~ cudaSetDevice(0);
	//~ if(d_resultsFinal)
	//~ {
		//~ cudaFree(d_resultsFinal);
	//~ }
	//~ CUERR
}

/**
 * Allocate pinned host memory and device memory for batch data and results.
 */
template <typename HKMERr>
size_t CuClarkDB<HKMERr>::malloc(size_t _numReads,
						size_t _maxReads, size_t _maxReadsInContainers,
						std::vector<ITYPE>& _indexBatches,
						RESULTS* &_fullResults,	size_t _resultRowSize,
						RESULTS* &_finalResults, size_t _finalResultsRowSize,
						bool _isExtended,
						std::vector<uint32_t*>& _readsPointer,
						std::vector<CONTAINER*>& _readsInCon)
{
	m_sizeResultRow			= _resultRowSize*sizeof(RESULTS);
	m_sizeResultFinalRow	= _finalResultsRowSize*sizeof(RESULTS);
	
	// set parameters for query kernel
	size_t numWarps = 2;
	size_t warpSize = 32;
	m_threadsPerBlock_queryKernel = warpSize*numWarps;
	//~ size_t numTargets = (((m_numTargets-1) / threadsPerBlock) +1) *threadsPerBlock;	
	
	// hit counters in shared
	m_sharedSize_queryKernel = ((m_numTargets-1)/2+1)*2*sizeof(uint16_t);
	// containers in shared
	size_t sharedMemPerWarp = ((m_k+warpSize-2)/(sizeof(CONTAINER)*4)+1)*sizeof(CONTAINER);
	// result buffer in shared
	if (m_sizeResultRow > sharedMemPerWarp) sharedMemPerWarp = m_sizeResultRow;
	m_sharedSize_queryKernel += sharedMemPerWarp*numWarps;

	//~ std::cerr << "Required Shared Memory per block: \t" << m_sharedSize_queryKernel/1000.0 << " KB\n";
	
	size_t total = 0;
	
	size_t maxReadsPointer = _maxReads+1;
	size_t sizeReadsPointer			= maxReadsPointer*sizeof(uint32_t),
		   sizeReadsInContainers	= _maxReadsInContainers*sizeof(CONTAINER);
		   //~ sizeResultRow 			= _resultRowSize*sizeof(RESULTS);
	
	for (int i=0; i<m_numBatches; i++)
	{
		cudaMallocHost(&h_readsPointer[i], sizeReadsPointer);
		cudaMallocHost(&h_readsInContainers[i], sizeReadsInContainers);
		CUERR
	}
	_readsPointer = h_readsPointer;
	_readsInCon = h_readsInContainers;
	
	for (int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		
		// allocate space for reads on each device
		cudaMalloc(&d_readsPointer[i], sizeReadsPointer);
		CUMEMERR
		cudaMalloc(&d_readsInContainers[i], sizeReadsInContainers);
		CUMEMERR	

		// allocate space for each partitial result & merging
		for (int j=0; j<d_results[i].size(); j++)
			cudaMalloc(&d_results[i][j], m_sizeResultRow*_maxReads);
		//~ std::cerr << i << " allocated\n";
		CUMEMERR
	}
	
	total += sizeReadsPointer;
	total += sizeReadsInContainers;
	total += m_sizeResultRow*_maxReads*d_results[0].size();
	
	// allocate space to store full results on host
	if (m_dbParts > 1 || _isExtended)
	{
		//~ std::cerr << "extra allocation\n";
		cudaMallocHost(&_fullResults, m_sizeResultRow*_numReads);
		CUERR
		//~ std::cerr << "Full result size:\t" << m_sizeResultRow*_numReads/1000/1000.0 << " MB\n";
	}
	
	// allocate space for final result & async copy to host
	if (_finalResultsRowSize > 0)
	{
		//~ cudaSetDevice(0);
		cudaMallocHost(&_finalResults, m_sizeResultFinalRow*_numReads);
		CUERR
		
		for (int i=0; i<m_numDevices; i++)
		{
			cudaSetDevice(i);
			
			cudaMalloc	  (&d_resultsFinal[i],  m_sizeResultFinalRow*_maxReads);
			CUMEMERR
		}
		
		total += m_sizeResultFinalRow*_maxReads;
	}
	//~ std::cerr << "Final result size:\t" << m_sizeResultFinalRow*_numReads/1000/1000.0 << " MB\n";
	
	for (int i=0; i<m_numBatches; i++)
	{
		h_results[i] = _fullResults+_resultRowSize*_indexBatches[i];	
		h_resultsFinal[i] = _finalResults+_finalResultsRowSize*_indexBatches[i];
	}
	
	return total;
}

/** 
 * Synchronize all CUDA devices and check for errors.
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::sync()
{
	for (int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
		CUERR
	}
	return true;
}

/** 
 * Wait for GPUs to finish a batch.
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::waitForBatch(size_t batchId)
{
	cudaEventSynchronize(m_batchFinishedEvents[batchId]);
	return cudaGetLastError() == cudaSuccess;
}

/** 
 * Check if GPUs have already finished a batch.
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::checkBatch(size_t batchId)
{
	cudaEventQuery(m_batchFinishedEvents[batchId]);
	return cudaGetLastError() == cudaSuccess;
}

/** 
 * Read database from files to pinned host memory,
 * divide into parts according to CUDA device memory.
 * cf. hashTable_hh read
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::read (const char * _filename, size_t& _fileSize, size_t& _dbParts, const ITYPE& _modCollision, const bool& _isfastLoadingRequested)
{

	char * file_sze = (char*) calloc(strlen(_filename)+4,sizeof(char));
	char * file_key = (char*) calloc(strlen(_filename)+4,sizeof(char));
	char * file_lbl = (char*) calloc(strlen(_filename)+4,sizeof(char));

	sprintf(file_sze, "%s.sz", _filename);
	sprintf(file_key, "%s.ky", _filename);
	sprintf(file_lbl, "%s.lb", _filename);

	std::ifstream ifs_sze;
	std::ifstream ifs_key;
	std::ifstream ifs_lbl;

//~ #define BUFFERSIZE 64001

	//~ char buffer_sze[BUFFERSIZE];
	//~ ifs_sze.rdbuf()->pubsetbuf(buffer_sze,BUFFERSIZE);
	//~ char buffer_key[BUFFERSIZE];
	//~ ifs_key.rdbuf()->pubsetbuf(buffer_key,BUFFERSIZE);
	//~ char buffer_lbl[BUFFERSIZE];	
	//~ ifs_lbl.rdbuf()->pubsetbuf(buffer_lbl,BUFFERSIZE);
	
	ifs_sze.open(file_sze, std::ios::binary);
	ifs_key.open(file_key, std::ios::binary);
	ifs_lbl.open(file_lbl, std::ios::binary);
	
	if (!ifs_sze.is_open())
	{	std::cerr << "Failed to open " << file_sze << std::endl; return false;	}
	if (!ifs_key.is_open())
	{	std::cerr << "Failed to open " << file_key << std::endl; return false;   }
	if (!ifs_lbl.is_open())
	{	std::cerr << "Failed to open " << file_lbl << std::endl; return false;   }
	
	bool allCollision = _modCollision <= 1;
	
	/// read bucket sizes
	std::vector<uint8_t>	bucketSizes;
	bucketSizes.resize(HTSIZE);
	
	ifs_sze.read((char*) &bucketSizes[0], HTSIZE);

	uint64_t nbElements = 0;
	uint32_t nbNonZeroBuckets = 0;
	//~ size_t bucketSizeMax = 0;
	std::vector<uint8_t>	choice(HTSIZE,0);
	
	// choose buckets and count chosen elements
	for (uint32_t i = 0; i < HTSIZE; i++)
	{
		if (bucketSizes[i] > 0)
		{
			nbNonZeroBuckets++;
			// 2 = keep, 1 = skip, 0 = empty
			choice[i] = (allCollision || (nbNonZeroBuckets % _modCollision)== 0) ? 2: 1;
			if (choice[i] == 2)
				nbElements += bucketSizes[i];
				
			//~ if (bucketSizeMax < bucketSizes[i])
				//~ bucketSizeMax = bucketSizes[i];
		}
	}
	
	//~ std::cerr << "AVG bucket size: " << (float)nbElements/HTSIZE
			  //~ << ", MAX bucket size: " << (int)bucketSizeMax;
	
	// calculate total size
	_fileSize = HTSIZE;	
	// bucket pointers: 8bit -> 32bit per element
	_fileSize *= sizeof(uint32_t);
	// size of keys and labels
	size_t _fileSizeKeys   = nbElements * sizeof(HKMERr);
	size_t _fileSizeLabels = nbElements * sizeof(ILBL);
	// total database size
	_fileSize = _fileSize + _fileSizeKeys + _fileSizeLabels;
	std::cerr << "Total DB size in RAM:\t" << _fileSize/1000000/1000.0 << " GB\n";
	
	/// divide into parts
	std::cerr << "Total device memory:\t" << m_memSizes.back()/1000000/1000.0 << " GB (" << m_numDevices*RESERVED/1000000 << " MB reserved)\n";
	
	size_t minParts = (nbElements/(uint32_t)-1)+1;
	//~ std::cerr << ", Need at least " << minParts << " part(s).\n";

	m_cyclesPerDevice = _fileSize / m_memSizes.back() + 1;
	m_dbParts = m_cyclesPerDevice*m_numDevices*m_dbPartsPerDevice;
	
	// adjust number of parts to prevent overflow
	if (m_dbParts < minParts)
	{
		std::cerr << "Overflow prevented.\n";
		m_cyclesPerDevice = (minParts-1)/(m_numDevices*m_dbPartsPerDevice)+1;
		m_dbParts = m_cyclesPerDevice*m_numDevices*m_dbPartsPerDevice;
	}
	
	m_cyclesToDo = m_cyclesPerDevice;
	_dbParts = m_cyclesPerDevice;
	
	// no space for merging needed if db fits on one device
	if (m_dbParts == 1)
		d_results[0].resize(1);
	
	std::cerr << "Requiring " << m_cyclesPerDevice << " loop(s).\n";

	// point to db parts in bucket list
	//~ std::vector<uint32_t> m_partPointer(m_dbParts+1, 0);		// named Pos in hashTable_hh
	m_partPointer.resize(m_dbParts+1);
	m_partPointer[0] = 0;
	for(int i = 0; i < m_dbParts ; i++)
	{
		//~ m_partPointer[i] = ((HTSIZE-1)/m_dbParts +1) * i;		
		//~ m_partPointer[i] = HTSIZE/m_dbParts * i;
		m_partPointer[1+i] = HTSIZE/m_cyclesPerDevice * (i/(m_numDevices*m_dbPartsPerDevice) + (float)m_memSizes[i%(m_numDevices*m_dbPartsPerDevice)] / m_memSizes[m_numDevices*m_dbPartsPerDevice-1]);
	}	
	m_partPointer[m_dbParts] = HTSIZE;
	
	//~ for(int i = 0; i <= m_dbParts ; i++)
	//~ {
		//~ std::cerr << "m_partPointer " << m_partPointer[i] << " \n";
	//~ }
	
	/// calculate pointers to buckets
	
	// pointers for each part
	h_bucketPointers.resize(m_dbParts);
	h_keys.resize(m_dbParts);
	h_labels.resize(m_dbParts);
	
	// sizes for each part
	m_partSize.resize(m_dbParts);
	m_partSizeKeys.resize(m_dbParts);
	m_partSizeLabels.resize(m_dbParts);
	
	/*	single part version
	std::vector<uint32_t>	h_bucketPointers;
	h_bucketPointers.resize(HTSIZE+1);
	h_bucketPointers[0] = 0;
	
	// partial_sum
	std::vector<uint8_t>::const_iterator first = bucketSizes.begin();
	std::vector<uint8_t>::const_iterator last  = bucketSizes.end();
	std::vector<uint32_t>::iterator result = h_bucketPointers.begin()+1;

	if (first!=last) {
	  uint32_t val = *first;
	  *result = val;
	  while (++first!=last) {
		val = val + *first;
		*++result = val;
	  }
	  ++result;
	}
	
	//~ std::copy(h_bucketPointers.begin(), h_bucketPointers.begin()+5,
				//~ std::ostream_iterator<unsigned>(std::cerr, " "));	          		
	//~ std::cerr << '\n';
	*/
		
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int i=0; i<m_dbParts; i++)
	{
		size_t numBuckets = m_partPointer[i+1]-m_partPointer[i];
		m_partSize[i] = (numBuckets + 1) * sizeof(uint32_t);
		cudaMallocHost(&h_bucketPointers[i], m_partSize[i]);
		CUERR
		
		h_bucketPointers[i][0] = 0;
		
		// partial_sum
		size_t firstIndex = m_partPointer[i];
		size_t lastIndex = m_partPointer[i+1];
		uint32_t* result = h_bucketPointers[i]+1;

		if (firstIndex < lastIndex)
		{
		  uint32_t val = 0;
		  if(choice[firstIndex] == 2)
		  {
			val = bucketSizes[firstIndex];
		  }
		  *result = val;
		  while (++firstIndex < lastIndex)
		  {
			if(choice[firstIndex] == 2)
			{
				val = val + bucketSizes[firstIndex];
			}
			if (val < *result)
			{
				std::cerr << "Bucket pointer overflow. Abort.\n";
			}
			*++result = val;
		  }
		  ++result;
		}
		
		m_partSizeKeys[i]   = h_bucketPointers[i][numBuckets] * sizeof(HKMERr);
		m_partSizeLabels[i] = h_bucketPointers[i][numBuckets] * sizeof(ILBL);

		//~ std::cerr << i+1 << "/" << (int)m_dbParts << " database: " << h_bucketPointers[i][numBuckets] << " elements, "
				  //~ << (m_partSize[i]+m_partSizeKeys[i]+m_partSizeLabels[i])/1000/1000.0 << " MB\n";
		
		//~ std::cerr << "Pointers: " << m_partSize[i]/1000/1000.0
				  //~ << " MB\t Keys: " << m_partSizeKeys[i]/1000/1000.0
				  //~ << " MB\t Labels: " << m_partSizeLabels[i]/1000/1000.0
				  //~ << " MB\n";
	}
	
	// point to db parts in keys and labels
	m_partPointerKeys.resize(m_dbParts+1);
	m_partPointerKeys[0] = 0;
	
	std::vector<size_t>	max_partSize(m_numDevices*m_dbPartsPerDevice,0),
						max_partSizeKeys(m_numDevices*m_dbPartsPerDevice,0),
						max_partSizeLabels(m_numDevices*m_dbPartsPerDevice,0);	
	
	for (int i=0; i<m_dbParts; i++)
	{
		size_t numBuckets = m_partPointer[i+1]-m_partPointer[i];
		m_partPointerKeys[i+1] = m_partPointerKeys[i] + h_bucketPointers[i][numBuckets];
		
		int device = (i/m_dbPartsPerDevice) % m_numDevices;
		//~ std::cerr << "max size - Device " << device << " Index " << i << std::endl;
		if (max_partSize[device] 		< m_partSize[i]) 		max_partSize[device] 		= m_partSize[i];
		if (max_partSizeKeys[device]  	< m_partSizeKeys[i]) 	max_partSizeKeys[device] 	= m_partSizeKeys[i];
		if (max_partSizeLabels[device]  < m_partSizeLabels[i]) 	max_partSizeLabels[device] 	= m_partSizeLabels[i];
	}
	//~ for (int i=0; i<m_dbParts; i++)
	//~ {
		//~ std::cerr << max_partSize[i] << " " << max_partSizeKeys[i] << " " << max_partSizeLabels[i] << "\n";
	//~ }
	
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	{	
		/// read bucket contents
		#ifdef _OPENMP
		#pragma omp single nowait
		#endif
		{
			for (int i=0; i<m_dbParts; i++)
			{				
				cudaHostAlloc(&h_keys[i], m_partSizeKeys[i], 0);
				//~ cudaHostAlloc(&h_labels[i], m_partSizeLabels[i], 0);
				CUERR
				
				//~ std::cerr << "Reading...\n";
				if (_modCollision <= 1)
				{	// read everything
					ifs_key.read((char*) &h_keys[i][0], m_partSizeKeys[i]);
					//~ ifs_lbl.read((char*) &h_keys[i][0], m_partSizeLabels[i]);
				}
				else
				{	// read if choice == 2
					uint8_t bucketSize;
					size_t ignoreSize = 0;
					//~ char keybuffer[10241];
					//~ ifs_key.rdbuf()->pubsetbuf(keybuffer,10241);
					uint32_t storeIndexKeys = 0;
					// ~ uint32_t storeIndexLabels = 0;
					for (uint32_t j = m_partPointer[i]; j < m_partPointer[i+1]; j++)
					{
						bucketSize = bucketSizes[j];
						switch(choice[j])
						{
							case 1:
								ignoreSize += bucketSize;
								// skip
								//~ ifs_key.ignore(bucketSize*sizeof(HKMERr));
								//~ ifs_lbl.ignore(bucketSize*sizeof(ILBL));
								break;
							case 2:
								//skip
								ifs_key.ignore(ignoreSize*sizeof(HKMERr));
								//~ ifs_lbl.ignore(ignoreSize*sizeof(ILBL));
								ignoreSize = 0;
								// read
								ifs_key.read((char*) &h_keys[i][storeIndexKeys], bucketSize*sizeof(HKMERr));
								storeIndexKeys += bucketSize;
								// ~ ifs_lbl.read((char*) &h_labels[storeIndexLabels], bucketSize*sizeof(ILBL));
								// ~ storeIndexLabels += bucketSize;
								break;
						}
					}
				}
			}
		}
		#ifdef _OPENMP
		#pragma omp single
		#endif
		{
			for (int i=0; i<m_dbParts; i++)
			{
				// ~ cudaHostAlloc(&h_keys[i], m_partSizeKeys[i], 0);
				cudaHostAlloc(&h_labels[i], m_partSizeLabels[i], 0);
				CUERR
			
				//~ std::cerr << "Reading...\n";
				if (_modCollision <= 1)
				{	// read everything
					//~ ifs_key.read((char*) &h_keys[i][0], m_partSizeKeys[i]);
					ifs_lbl.read((char*) &h_labels[i][0], m_partSizeLabels[i]);
				}
				else
				{	// read if choice == 2
					uint8_t bucketSize;
					size_t ignoreSize = 0;				
					//~ char lblbuffer[10241];	
					//~ ifs_lbl.rdbuf()->pubsetbuf(lblbuffer,10241);
					//~ uint32_t storeIndexKeys = 0;
					uint32_t storeIndexLabels = 0;
					for (uint32_t j = m_partPointer[i]; j < m_partPointer[i+1]; j++)
					{
						bucketSize = bucketSizes[j];
						switch(choice[j])
						{
							case 1:
								ignoreSize += bucketSize;
								// skip
								//~ ifs_key.ignore(bucketSize*sizeof(HKMERr));
								//~ ifs_lbl.ignore(bucketSize*sizeof(ILBL));
								break;
							case 2:
								// skip
								//~ ifs_key.ignore(ignoreSize*sizeof(HKMERr));
								ifs_lbl.ignore(ignoreSize*sizeof(ILBL));
								ignoreSize = 0;
								// read
								//~ ifs_key.read((char*) &h_keys[i][storeIndexKeys], bucketSize*sizeof(HKMERr));
								//~ storeIndexKeys += bucketSize;
								ifs_lbl.read((char*) &h_labels[i][storeIndexLabels], bucketSize*sizeof(ILBL));
								storeIndexLabels += bucketSize;
								break;
						}
					}
				}
			}
		}

		#ifdef _OPENMP
		#pragma omp single
		#endif	
		for (int i=0; i<m_numDevices; i++)
		{
			cudaSetDevice(i);
			
			for(int j=0; j<m_dbPartsPerDevice; ++j)
			{
				int index = m_dbPartsPerDevice*i+j;
				//~ std::cerr << "alloc - Device " << i << " Index " << index << std::endl;

				// allocate device memory
				cudaMalloc(&d_bucketPointers[index], max_partSize[i]);
				CUERR
				cudaMalloc(&d_keys[index], max_partSizeKeys[i]);
				CUERR
				cudaMalloc(&d_labels[index], max_partSizeLabels[i]);
				CUERR
			}
		}
 	}
 	std::cerr << "DB loaded in RAM.\n";
     
	return true;
}

/**
 *  Swap to next database parts on all GPUs
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::swapDbParts ()
{
	if (m_cyclesToDo == 0)
	{
		// reset for possible next file
		m_cyclesToDo = m_cyclesPerDevice;
		return false;
	}
	
	//~ if (m_cyclesToDo < m_cyclesPerDevice)
		//~ std::cerr << "Swapping DB part. ";
	
	int offset = m_numDevices*m_dbPartsPerDevice*(m_cyclesPerDevice - m_cyclesToDo);
	//~ std::cerr << "Offset: " << offset << "\n";
	
	//~ #ifdef _OPENMP
	//~ #pragma omp parallel for
	//~ #endif
	for (int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		//~ cudaDeviceSynchronize();
		//~ CUERR
		for(int j=0; j<m_dbPartsPerDevice; ++j)
		{
			int index = m_dbPartsPerDevice*i+j;
			//~ std::cerr << "Swap - Device " << i << " Index " << index << " Offset " << offset << std::endl;
			
			// copy database to part device
			cudaMemcpyAsync(d_bucketPointers[index], &h_bucketPointers[index+offset][0], m_partSize[index+offset], cudaMemcpyHostToDevice, 0);
			//~ CUERR
			cudaMemcpyAsync(d_keys[index],   &h_keys  [index+offset][0], m_partSizeKeys[index+offset],   cudaMemcpyHostToDevice, 0);
			//~ CUERR
			cudaMemcpyAsync(d_labels[index], &h_labels[index+offset][0], m_partSizeLabels[index+offset], cudaMemcpyHostToDevice, 0);
			//~ CUERR
		}
		
		//~ cudaDeviceSynchronize();
		//~ CUERR
	}
	//~ if (m_cyclesToDo < m_cyclesPerDevice)
		//~ std::cerr << "DB parts to do: " << m_cyclesToDo << std::endl;
	
	m_cyclesToDo--;
	return true;
}

/**
 *  Set all parameters for a batch
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::readyBatch (const size_t _batchId, const size_t _numReads, const size_t _containerCount)
{
	m_numReads[_batchId] = _numReads;
	m_sizeReadsPointer[_batchId] = (_numReads+1)*sizeof(uint32_t);
	m_sizeReadsInContainers[_batchId] = _containerCount*sizeof(CONTAINER);

	return true;
}

/**
 * 	Schedule a prepared batch for database query on GPUs,
 *  Schedule results merging if needed,
 * 	Schedule calculation of final results at the end.
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::queryBatch (const size_t _batchId, const bool _isExtended, const bool _isFollowup)
{
	size_t d_pitch = m_sizeResultRow;
	
	cudaStream_t stream = 0;
    //~ cudaStreamCreate(&stream);
    
	for (int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		cudaMemcpyAsync(d_readsPointer[i], h_readsPointer[_batchId], m_sizeReadsPointer[_batchId], cudaMemcpyHostToDevice, stream);
		//~ CUERR
		cudaMemcpyAsync(d_readsInContainers[i], h_readsInContainers[_batchId], m_sizeReadsInContainers[_batchId], cudaMemcpyHostToDevice, stream);
		//~ CUERR
	}
	
	// process one read per block
	size_t numBlocks = m_numReads[_batchId];
	
	int dbPartOffset = m_numDevices*m_dbPartsPerDevice*(m_cyclesPerDevice - m_cyclesToDo -1);
	//~ std::cerr << "DB Offset: " << dbPartOffset << "\n";
	
	for (int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		
		for(int j=0; j<m_dbPartsPerDevice; ++j)
		{
			int index = m_dbPartsPerDevice*i+j;
			//~ if(_batchId==0) std::cerr << "query - Device " << i << " Index " << index << " Offset " << dbPartOffset << std::endl;
			queryKernel<<<numBlocks, m_threadsPerBlock_queryKernel, m_sharedSize_queryKernel, stream>>>
							(m_k,
							d_readsPointer[i], d_readsInContainers[i],
							d_bucketPointers[index], d_keys[index], d_labels[index],
							m_partPointer[index+dbPartOffset], m_partPointer[index+dbPartOffset+1],
							d_results[i][j], d_pitch, m_numTargets);
			//~ cudaStreamSynchronize(stream);
			//~ CUERR
		}
	}
		
	size_t 	threadsPerBlock = 1024;		
	numBlocks = (m_numReads[_batchId]-1) / threadsPerBlock +1;
	
	// merge results from different db parts on same device
	for (int i=0; i<m_numDevices; i++)
	{
		cudaSetDevice(i);
		
		for(int j=1; j<m_dbPartsPerDevice; j<<=1)
		{
			for (int k=0; j+k<m_dbPartsPerDevice; k+=2*j)
			{
				//~ if(_batchId==0) std::cerr << "merge - Device " << i << std::endl;

				mergeKernel<<<numBlocks,threadsPerBlock,0,stream>>>(d_results[i][k], d_results[i][k+j], d_pitch, m_numReads[_batchId], d_results[i][m_dbPartsPerDevice]);
				//~ cudaStreamSynchronize(stream);
				//~ CUERR
				
				// swap pointers so that merged result is first
				RESULTS* dummy = d_results[i][k];
				d_results[i][k] = d_results[i][m_dbPartsPerDevice];
				d_results[i][m_dbPartsPerDevice] = dummy;
			}
		}
	}	
	
	// merge results from different devices
	for (int i=1; i<m_numDevices; i<<=1)
	{
		for (int j=0; j+i<m_numDevices; j+=2*i)
		{
			cudaSetDevice(j);
			// sync devices, async to host
			cudaMemcpyPeer(d_results[j][1], j, d_results[j+i][0], j+i, m_sizeResultRow*m_numReads[_batchId]);	
			// async to other streams
			//~ cudaMemcpyPeerAsync(d_results[2], 0, d_results[1], 1, m_sizeResultRow*m_numReads[_batchId], stream);	
			
			mergeKernel<<<numBlocks,threadsPerBlock,0,stream>>>(d_results[j][0], d_results[j][1], d_pitch, m_numReads[_batchId], d_results[j][2]);
			//~ cudaStreamSynchronize(stream);
			//~ CUERR
			
			// swap pointers so that merged result is first
			RESULTS* dummy = d_results[j][0];
			d_results[j][0] = d_results[j][2];
			d_results[j][2] = dummy;			
		}
	}
	
	cudaSetDevice(0);
	// merge with previous results
	if(_isFollowup)
	{
		//~ std::cerr << "Batch " << _batchId << " retrieve from host" << std::endl;
		
		cudaMemcpyAsync(d_results[0][1], h_results[_batchId], m_sizeResultRow*m_numReads[_batchId], cudaMemcpyHostToDevice, stream);
				
		mergeKernel<<<numBlocks,threadsPerBlock,0,stream>>>(d_results[0][0], d_results[0][1], d_pitch, m_numReads[_batchId], d_results[0][2]);
		//~ cudaStreamSynchronize(stream);
		//~ CUERR
		
		// swap pointers so that merged result is first
		RESULTS* dummy = d_results[0][0];
		d_results[0][0] = d_results[0][2];
		d_results[0][2] = dummy;		
	}
	
	// copy result to host
	if (m_cyclesToDo > 0 || _isExtended)
	{	// copy full results to host
		//~ std::cerr << "Batch " << _batchId << " copy to host" << std::endl;
				
		cudaMemcpyAsync(h_results[_batchId], d_results[0][0], m_sizeResultRow*m_numReads[_batchId], cudaMemcpyDeviceToHost, stream);
		//~ cudaStreamSynchronize(stream);
		//~ CUERR

		//~ std::cerr << "Result size: " << m_sizeResultRow*_numReads / 1000000.0 << " MB" << std::endl;
	}
	
	if (m_cyclesToDo == 0)
	{	// calculate final results and copy to host
		//~ std::cerr << "Batch " << _batchId << " calculate results" << std::endl;
		
		resultKernel<<<numBlocks,threadsPerBlock,0,stream>>>(d_results[0][0], d_pitch, m_numReads[_batchId], d_resultsFinal[0], m_sizeResultFinalRow);
		//~ cudaStreamSynchronize(stream);
		//~ CUERR
		
		cudaMemcpyAsync(h_resultsFinal[_batchId], d_resultsFinal[0], m_sizeResultFinalRow*m_numReads[_batchId], cudaMemcpyDeviceToHost, stream);
		//~ cudaStreamSynchronize(stream);
		//~ CUERR
		
		//~ std::cerr << "Result size: " << m_sizeResultFinalRow*m_numReads[_batchId] / 1000000.0 << " MB" << std::endl;

		cudaEventRecord(m_batchFinishedEvents[_batchId], stream);
		return true;
	}
	
	//~ cudaStreamSynchronize(stream);
	//~ CUERR	
	//~ cudaDeviceSynchronize();
	//~ CUERR
	return false;
}

/**
 * Get the two targets with highest scores and sum of all targets.
 */
template <typename HKMERr>
bool CuClarkDB<HKMERr>::getFinalResult (const size_t _batchId, RESULTS* _resutlsFinal)
{
	size_t d_pitch = m_sizeResultRow;
	
	int device = _batchId % m_numDevices;
	cudaSetDevice(device);
	
	cudaMemcpyAsync(d_results[device][0], h_results[_batchId], m_sizeResultRow*m_numReads[_batchId], cudaMemcpyHostToDevice, 0);
	
	size_t threadsPerBlock = 1024;
	size_t numBlocks = (m_numReads[_batchId]-1) / threadsPerBlock +1;
	//~ std::cerr << "# reads: " << m_numReads[_batchId] << " # blocks: " << numBlocks << "\n";
	resultKernel<<<numBlocks,threadsPerBlock,0,0>>>(d_results[device][0], d_pitch, m_numReads[_batchId], d_resultsFinal[device], m_sizeResultFinalRow);
	
	cudaMemcpyAsync(_resutlsFinal, d_resultsFinal[device], m_sizeResultFinalRow*m_numReads[_batchId], cudaMemcpyDeviceToHost, 0);
	
	cudaEventRecord(m_batchFinishedEvents[_batchId], 0);
	
	return true;
}


/**
 * Queries a batch against a database part.
 * 
 * Processes one read per block:
 * Loads read data into shared memory,
 * contructs one kmer per thread and queries it,
 * scores in shared memory,
 * continues until the read is completed.
 * Stores non zero scores in global memory.
 */
template <typename HKMERr>
__global__ void queryKernel (uint8_t k,
			uint32_t* readsPointer, CONTAINER* readsInContainers,
			uint32_t* bucketPointers, HKMERr* keys, ILBL* labels,
			uint32_t dbPartStart, uint32_t dbPartEnd,
			RESULTS* results, size_t pitch, size_t numTargets)
{
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int wid = threadIdx.x/warpSize;
	int wlane = threadIdx.x % warpSize;
	int numWarps = blockDim.x/warpSize;
	size_t nucsPerCon = sizeof(CONTAINER)*4;
	
	// get row (=bid) to store results
	RESULTS* resultRow = (RESULTS*) ((char*)results + bid*pitch);
	//~ uint32_t* resultRow = &results[bid*numTargets];
	
	// shared memory for scoring
	extern __shared__ uint32_t targetHits[];
	uint16_t* targetHits16 = (uint16_t*) &targetHits[0];
	RESULTS* sharedResultRow = (RESULTS*) &targetHits16[((numTargets-1)/2+1)*2];
	CONTAINER* sharedContainers = (CONTAINER*) &sharedResultRow[0];
	
	// set all target counters to zero
	for (int i=tid; i<numTargets; i += blockDim.x)
	{
		targetHits16[i] = 0;
	}
	__syncthreads();	// for blockDim > warp size
	
	uint64_t kmer;
	ILBL target;
	uint64_t cutoff = (uint64_t)-1 >> (64 - 2*k);
	
	CONTAINER partLength, offset, tmp;	
	uint32_t partPointer = readsPointer[bid];
	uint32_t readEnd = readsPointer[bid+1];
	//~ uint32_t partBegin;
	uint32_t partIterator;
	short numKmer;
	
	uint32_t firstContainer;
	uint32_t containerPerWarp = (k+warpSize-2)/nucsPerCon+1;
	//~ if(tid==0) printf("Block %d: MaxCon %u\n",bid,maxContainers);
	
	while(partPointer < readEnd)
	{
		partLength = readsInContainers[partPointer];
		//~ partBegin = ++partPointer + tid / nucsPerCon;
		firstContainer = ++partPointer;
		partPointer += (partLength-1) / nucsPerCon +1;
		// number of kmer for wlane == 0 to check if warp has work to do
		numKmer = partLength - k + 1 - wid*warpSize;
		while (numKmer > 0)
		{
			// load containers for own warp
			if(wlane < containerPerWarp)
			{
				int readIndex = firstContainer+wid*warpSize/nucsPerCon+wlane;
				if (readIndex < partPointer)
					sharedContainers[wid*containerPerWarp + wlane] 
						= readsInContainers[readIndex];
			}
			// check if thread has work to do
			if(numKmer-wlane <= 0) break;
			partIterator = wid*containerPerWarp + wlane / nucsPerCon;	// shared version
			//~ partIterator = firstContainer + tid / nucsPerCon;	// global version
			
			kmer = 0;
			// read full containers
			for (int i=0; i<(k+tid%nucsPerCon)/ nucsPerCon; ++i)
			{
				kmer <<= 2*nucsPerCon;
				kmer |= sharedContainers[partIterator++];	// shared version
				//~ kmer |= readsInContainers[partIterator++];	// global version
			}
			// read 'offset' additional nucs
			offset = (k+tid) % nucsPerCon;
			if (offset != 0)
			{				
				kmer <<= 2*offset;
				tmp = sharedContainers[partIterator];	// shared version
				//~ tmp = readsInContainers[partIterator];	// global version
				tmp >>= 2*(nucsPerCon-offset);
				kmer |= tmp;
			//~ if(tid==warpSize-1) printf("Block %d: Iterator %u\n",bid,partIterator);
			}
			// cut off overhang
			kmer &= cutoff;
			
			// print for debugging
			/*
			printf("Block %d, Thread %2d: %" PRIu64 "\n",bid,tid,kmer);
			uint64_t x = 3;
			char kmer_string[28];
			kmer_string[27] = '\0';
			for(int j=0; j<k; ++j)
			{
				switch ( (kmer & x) >> (2*j) )
				{
					case(0): kmer_string[k-j-1] = 'T'; break;
					case(1): kmer_string[k-j-1] = 'G'; break;
					case(2): kmer_string[k-j-1] = 'C'; break;
					case(3): kmer_string[k-j-1] = 'A'; break;
				}
				x <<= 2;
			}
			printf("Block %d, Thread %2d: %s\n",bid,tid,kmer_string);
			*/
			
			if (queryElement(k, kmer, bucketPointers, keys, labels, dbPartStart, dbPartEnd, target))
			{
				// 32bit targetHits version
				//~ atomicAdd(&targetHits[target], 1);
				//~ printf("Target: %d\n",target);
				
				// 16bit targetHits version
				ILBL target32 = target / 2;
				//~ uint32_t value = target % 2 ? 1<<16 : 1;
				uint32_t value = 1 <<((target % 2)*16);
				atomicAdd(&targetHits[target32],value);
				//~ printf("Block: %d, Target %d, counter %d\n",bid, target,targetHits16[target]);
				//~ printf("Target: %d, target32: %d, value %d, counter: %d\n",target,target32,value,targetHits16[target]);
			}
			// advance to next section
			numKmer -= blockDim.x;
			firstContainer += blockDim.x/nucsPerCon;
		}				
	}
	__syncthreads();	// for blockDim > warp size
	
	// store nonzeros with index in global
	// cf. http://www.davidespataro.it/cuda-stream-compaction-efficient-implementation/
	//~ if (wid == 0)
	{
		int pred, t_m, b, t_u, total=0;
		int j;
		size_t numTargetsPerWarp = (numTargets-1) / numWarps +1;
		numTargetsPerWarp = ((numTargetsPerWarp-1) / warpSize +1) * warpSize;
		//~ if(wlane==0) printf("numTargetsPerWarp: %d\n",numTargetsPerWarp);
		//~ for (int i=tid; i<numTargets; i += warpSize)
		for (int i=wlane+numTargetsPerWarp*wid; (i<numTargetsPerWarp*(wid+1)) && (i <numTargets) ; i += warpSize)
		{
			pred = targetHits16[i] > 0 ? 1 : 0;
			t_m = INT_MAX >> warpSize-wlane-1;	// set bits < tid
			b = __ballot(pred) & t_m;			// get pred bits < tid
			t_u = __popc(b);					// get sum of bits = # pred < tid
			j = 2*(total+t_u);

			if (pred)
			{
				if (j+2 < pitch)
				{
					sharedResultRow[j+1+wid*(pitch/sizeof(RESULTS))] = i;
					sharedResultRow[j+2+wid*(pitch/sizeof(RESULTS))] = targetHits16[i];
					//~ printf("Block: %i, Target: %i, Hits: %i\n",bid,i,targetHits[i]);
				}
				else
				{
					printf("Too many different tagets hit by a sequence. Results will be corrupted.\n");
				}
			}
			total += t_u+pred;
			if (i == numTargetsPerWarp*(wid+1)-1 || i == numTargets-1)
			{
				//~ printf("Block: %i, Total: %i\n",bid,total);
				sharedResultRow[0+wid*(pitch/sizeof(RESULTS))] = total;
			}
			total = __shfl(total, warpSize-1);	// get total from last lane
		}
		

	}
	__syncthreads();	// for blockDim > warp size	
	
	int subtotal = 0;
	int warpTotal = sharedResultRow[wid*(pitch/sizeof(RESULTS))];
	for(int j=0; j<wid; j++)
	{
		subtotal += sharedResultRow[j*(pitch/sizeof(RESULTS))];	
	}
	for (int i=wlane; i<warpTotal*2; i+= warpSize)
			resultRow[i+subtotal*2+1] = sharedResultRow[i+wid*(pitch/sizeof(RESULTS))+1];
	if(wid==numWarps-1 && wlane == 0)
	{
		resultRow[0] = subtotal+warpTotal;
		//~ printf("Block: %d, Targets hit: %d\n",bid,subtotal+warpTotal);
	}
}

/**
 *  Query k-mer against a database part.
 *  Analog to hashTable_hh find, for canonical kmer
 */
template <typename HKMERr>
__device__ bool queryElement (const uint8_t& k, const uint64_t& _ikmer,
		uint32_t* d_bucketPointers, HKMERr* d_keys, ILBL* d_labels,
		uint32_t dbPartStart, uint32_t dbPartEnd,
		ILBL& _returnLabel)
{
	// getting reverse kmer
	size_t _ikmerR = _ikmer;
	// The following 6 lines come from Jellyfish source code
	_ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
	_ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
	_ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
	_ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
	_ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
	_ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (k << 1));
	
	// getting canonical kmer
	size_t _ikmerC = _ikmer < _ikmerR ? _ikmer : _ikmerR;
	
	size_t quotient = _ikmerC / HTSIZE;
	size_t remainder = _ikmerC - quotient * HTSIZE;
	
	// check for correct dbPart
	if (remainder < dbPartStart || remainder >= dbPartEnd)
		return false;
	remainder -= dbPartStart;
	
	//~ if(dbPart==1)
		//~ printf("Part: %d, Kmer: %" PRIu64 ", Remainder: %" PRIu64 "\n", dbPart, _ikmerC, remainder);

	size_t bucketBegin = d_bucketPointers[remainder];
	size_t bucketEnd = d_bucketPointers[remainder+1];
	//~ if(dbPart==1)
		//~ printf("Bucket size: %2d, Remainder: %" PRIu64 "\n", bucketEnd-bucketBegin, remainder-HTSIZE/2);

	if(	bucketEnd-bucketBegin > 0)
	{	// bucket not empty	
		size_t i = bucketBegin;
		HKMERr key = d_keys[i];
		
		if( key > quotient || d_keys[bucketEnd-1] < quotient)
		{	// quotient not in range
			return false;
		}
		
		while(key <= quotient)
		{
			if(key == quotient)
			{	// key found
				_returnLabel = d_labels[i];
				//printf("Part: %d, Label: %4d, Remainder: %8d\n", dbPart, _returnLabel, remainder);
				return true;
			}
			key = d_keys[++i];
		}
		// key not in list
		return false;
	} 

	// bucket empty
	return false;
}

/**
 * Merge two results into one.
 * 1 thread handles results of 1 read.
 * cf. https://nvlabs.github.io/moderngpu/merge.html
 */
__global__ void mergeKernel (RESULTS* resultA, RESULTS* resultB, size_t pitch, size_t numReads, RESULTS* results)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (tid < numReads)
	{
		RESULTS* resultRowA = (RESULTS*) ((char*)resultA + tid*pitch);
		RESULTS* resultRowB = (RESULTS*) ((char*)resultB + tid*pitch);
		RESULTS* resultsRow = (RESULTS*) ((char*)results + tid*pitch);
		
		// get number of results
		RESULTS countA = resultRowA[0];
		RESULTS countB = resultRowB[0];
		int count = countA + countB;
		//~ printf("Thread: %d, CountA: %d, CountB: %d\n",tid,countA,countB);
		if (count > 4*MAXHITS+4)
		{
			printf("Count overflow\n");
			count = 0;
		}
		
		int indexA = 1;
		int indexB = 1;
		
		// get target ids of first results
		RESULTS targetA = resultRowA[indexA];
		RESULTS targetB = resultRowB[indexB];
		
		//~ for(int i=0; i<countA; ++i)
		//~ {
			//~ printf("Thread: %d, targetA: %d\n",tid,resultRowA[2*i+1]);
		//~ }
		//~ for(int i=0; i<countB; ++i)
		//~ {
			//~ printf("Thread: %d, targetB: %d\n",tid,resultRowB[2*i+1]);
		//~ }
		
		// merge results
		for(int i=0; i<count; ++i)
		{
			//~ printf("Thread: %d, targetA: %d, targetB: %d\n",tid,targetA,targetB);
			
			// find next result
			uint8_t choice;
			// choice 0 -> A
			// choice 1 -> B
			// choice 2 -> A=B, merge
			
			if (indexB > countB*2)		// B done?
				choice = 0;	
			else if (indexA > countA*2)	// A done?
				choice = 1;
			else if (targetA < targetB) // A smaller?
				choice = 0;
			else if (targetA > targetB)	// B smaller?
				choice = 1;
			else						// same value
				choice = 2;
			
			// put target id
			resultsRow[2*i+1] = choice ? targetB : targetA;
			// put target counter
			resultsRow[2*i+2] = (choice ? resultRowB[++indexB] : resultRowA[++indexA])
							+ (choice == 2 ? resultRowA[++indexA] : 0);
			
			// get next target ids
			switch(choice)
			{	
				case 0:
					targetA = resultRowA[++indexA];
					break;
				case 1:
					targetB = resultRowB[++indexB];
					break;
				case 2:
					targetA = resultRowA[++indexA];
					targetB = resultRowB[++indexB];
					// one less output because of sum
					count--;
					break;				
			}
		}
		// put number of results
		resultsRow[0] = count;

		//~ for(int i=0; i<count; ++i)
			//~ printf("Thread: %d, target: %d, hits: %d\n",tid,resultsRow[2*i+1],resultsRow[2*i+2]);
	}
}

/**
 * Find best, second best and sum of scores for each read.
 * 1 thread handles results of 1 read.
 */
__global__ void resultKernel (RESULTS* scores, size_t spitch, size_t numReads, RESULTS* results, size_t rpitch)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (tid < numReads)
	{
		RESULTS* scoresRow = (RESULTS*) ((char*)scores + tid*spitch);
		RESULTS* resultsRow = (RESULTS*) ((char*)results + tid*rpitch);
		
		RESULTS best = 0, s_best = 0;
		RESULTS indexBest = 0, index_sBest = 0;
		RESULTS sumN = 0;
		
		RESULTS count = scoresRow[0];
		//~ printf("Read: %d, count: %d\n", tid, count);
		
		RESULTS targetScore;
		
		for(int i=0; i<count; ++i)
		{
			targetScore = scoresRow[2*i+2];
			
			// new best, update best and second best
			if (targetScore > best)
			{
				s_best = best;
				index_sBest = indexBest;
				best = targetScore;
				indexBest = scoresRow[2*i+1] + 1;					
			}
			// new second best, update
			else if (targetScore > s_best)
			{
				s_best = targetScore;
				index_sBest = scoresRow[2*i+1] + 1;
			}
			sumN += targetScore;
		}
		
		resultsRow[0] = sumN;
		resultsRow[1] = indexBest;
		resultsRow[2] = best;
		resultsRow[3] = index_sBest;
		resultsRow[4] = s_best;
		
		//~ printf("Thread: %d, sum: %d, best: %d, sbest: %d\n",tid,sumN,indexBest,index_sBest);
	}
}

// instantiations
template class CuClarkDB<uint16_t>;
template class CuClarkDB<uint32_t>;
template class CuClarkDB<uint64_t>;
