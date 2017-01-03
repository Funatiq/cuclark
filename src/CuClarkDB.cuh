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
 
 #ifndef CUCLARKDB_
#define CUCLARKDB_

#include <vector>
#include <cuda_runtime.h>
#include "./dataType.hh"

template <typename HKMERr>
class CuClarkDB
{
	private:
		uint8_t		m_k;			// kmer size
		size_t		m_numTargets;	// targets in db	
		size_t		m_numBatches;
		int			m_numDevices;
		
		std::vector<size_t>	m_memSizes;
		
		// db
		int			m_dbParts;
		int			m_dbPartsPerDevice;
		int			m_cyclesPerDevice;
		int			m_cyclesToDo;
		
		// host pointers
		std::vector< uint32_t* >	h_bucketPointers;
		std::vector< HKMERr* >		h_keys;
		std::vector< ILBL* >		h_labels;
		
		std::vector<size_t>		m_partSize;
		std::vector<size_t>		m_partSizeKeys;
		std::vector<size_t>		m_partSizeLabels;
		std::vector<uint32_t>	m_partPointer;
		std::vector<uint64_t>	m_partPointerKeys;
			
		std::vector<uint32_t*>	h_readsPointer;
		std::vector<CONTAINER*>	h_readsInContainers;
		
		std::vector<size_t>		m_numReads;
		std::vector<size_t>		m_sizeReadsPointer;
		std::vector<size_t>		m_sizeReadsInContainers;
		
		std::vector<RESULTS*>	h_results;
		std::vector<RESULTS*>	h_resultsFinal;
		
		size_t					m_sizeResultRow;
		size_t					m_sizeResultFinalRow;
		
		// device pointers
		std::vector<uint32_t*>	d_bucketPointers;
		std::vector<HKMERr*>	d_keys;
		std::vector<ILBL*>		d_labels;
		
		std::vector<uint32_t*>	d_readsPointer;
		std::vector<CONTAINER*>	d_readsInContainers;
		
		std::vector<std::vector<RESULTS*> >	d_results;	
		RESULTS*				d_resultsFinal;
		
		// events
		std::vector<cudaEvent_t> m_batchFinishedEvents;

		// kernel parameters
		size_t 		m_threadsPerBlock_queryKernel;
		size_t		m_sharedSize_queryKernel;
	
	public:
		CuClarkDB();
		CuClarkDB(	const size_t _numDevices, 
					const uint8_t _k,
					const size_t _numBatches,
					const size_t _numTargets
					);
		
		~CuClarkDB();
		
		void free();
		
		size_t malloc(	size_t _numReads,
						size_t _maxReads,
						size_t _maxReadsInContainers,
						std::vector<ITYPE>& _indexBatches,
						RESULTS* &_fullResults,
						size_t _resultRowSize,
						RESULTS* &_finalResults,
						size_t _finalResultsRowSize,
						bool _isExtended,
						std::vector<uint32_t*>& _readsPointer,
						std::vector<CONTAINER*>& _readsInCon
						);
						
		bool sync();
		
		bool waitForBatch(size_t i);
		
		bool checkBatch(size_t i);
		
		bool read(	const char * 	_filename,
					size_t& 		_fileSize,
					size_t& 		_dbParts,
					const ITYPE& 	_modCollision = 1,
					const bool& 	_isfastLoadingRequested = false
					);
				
		bool swapDbParts();
				
		bool readyBatch(const size_t _batchId,
						const size_t _numReads,
						const size_t _containerCount
						);
						
		bool queryBatch(const size_t _batchId,
						const bool _isExtended,
						const bool _isFollowup=false
						);
						
		bool getFinalResult(const size_t _batchId,
							RESULTS* _finalResult
							);
};

#endif
