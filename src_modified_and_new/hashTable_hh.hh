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
 * Changed hash table operations to canonical k-mers.
 */

#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include<vector>
#include<stdint.h>
#include<string>
#include "./dataType.hh"
#include "stdint.h"

template <typename HKMERr, typename ELMTr> struct htCell
{
	HKMERr CKey;
	ELMTr CElement;

	htCell(): CElement()
	{}

	htCell(const size_t& _quotient, const ILBL& _label, const size_t& _count): CElement() 
	{
		CKey = _quotient;
		CElement.Set(_label,_count);
	}

	htCell<HKMERr, ELMTr>& operator=(const size_t& _q)
	{
		CKey = _q;
		return *this;
	}

	bool operator<(const htCell<HKMERr,ELMTr>& a) const
	{
		return (CKey < a.CKey);
	}

};

template <typename HKMERr, typename ELMTr>
class hTable
{
	private:
		std::vector< sVector< htCell<HKMERr,ELMTr> > >   	m_table;
		size_t							m_load;
		size_t							m_it_x;
		size_t							m_it_y;
		uint8_t							m_k;

	public:
		hTable();
		hTable(const uint8_t _k);

		~hTable();

		size_t Load() const 
		{	return m_load;	}

		void sortall(const size_t& _iteratorPos = 0);

		bool insert(const uint64_t& 			_kmer,
				const htCell<HKMERr,ELMTr>& 	_cell
			   );

		bool insert(const uint64_t& 			_kmer, 
				const ILBL& 			_label
			   );

		bool insert(const uint64_t& 			_kmer, 
				const ILBL& 			_label,
				const size_t& 			_count
			   );

		bool find(const uint64_t& 			_kmer, 
				size_t& 			_xElement, 
				size_t& 			_yElement, 
				ILBL& 				_label, 
				IOCCR& 				_mult, 
				ICount& _count
			 ) const;

		bool find(const IKMER& 				_ikmer,
				const size_t& 			_reminderI, 
				const size_t& 			_quotientI, 
				ILBL& _label
			 ) const;

		bool find(const uint64_t&                       _ikmer,
				ILBL&                           _label
			 );

		bool find(const uint64_t& 			_ikmer, 
				const uint64_t& 		_ikmerR, 
				ILBL& 				_label
			 ) const;

		void updateElement(const size_t&	 	_xElement, 
				const size_t& 			_yElement, 
				const size_t& 			_count, 
				const bool& 			_lbl, 
				const bool& 			_isSameLbl = false
				); 

		void clear();

		uint64_t write(const char*	 		_fileht, 
				const size_t& 			iteratorPos, 
				const bool& 			_clearAfter = true
			      );

		bool read(const char * 				_filename, 
				size_t& 			_fileSize, 
				const size_t&	 		_nbCPU = 1,
				const ITYPE&                    _modCollision = 1,
				const bool& 			_isfastLoadingRequested = false
			 );

		bool resetIterator();
		bool nextIterator();
		void markElementAtIterator();
		void updateLabelAtIterator(const ILBL& _label);
		bool next();
		void elementIterator(ELMTr& _cElement) const;
		void elementIterator(uint64_t& 			_kmer, 
				ELMTr& 				_cElement
			) const;

};

#endif

/* @Author: Rachid Ounit
 * @Name: hash table class
 *
 */
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

using namespace std;

	template <typename HKMERr, typename ELMTr>
hTable<HKMERr, ELMTr>::hTable(): m_load(0), m_it_x(0), m_it_y(0), m_k(0)
{
	m_table.resize(HTSIZE);
}
	template <typename HKMERr, typename ELMTr>
hTable<HKMERr, ELMTr>::hTable(const uint8_t _k): m_load(0), m_it_x(0), m_it_y(0), m_k(_k)
{
	m_table.resize(HTSIZE);
}

	template <typename HKMERr, typename ELMTr>
hTable<HKMERr, ELMTr>::~hTable()
{
}

	template <typename HKMERr, typename ELMTr>
void hTable<HKMERr, ELMTr>::clear()
{
	m_load = 0;
	m_table.clear();
}

	template <typename HKMERr, typename ELMTr>
void hTable<HKMERr, ELMTr>::sortall(const size_t& _iteratorPos)
{
	size_t maxSize = 0;
	for(size_t t = 0; t < m_table.size() ; t++)
	{
		maxSize = !m_table[t].empty() ? (m_table[t].size() > maxSize? m_table[t].size(): maxSize): maxSize;
		if (!m_table[t].empty())
		{
			std::sort(m_table[t].begin()+_iteratorPos, m_table[t].end());
		}
	}
	cerr << "Hashtable sorting done: maximum number of collisions: " << maxSize -_iteratorPos << endl;
}

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::insert(const uint64_t& _kmer, const ILBL& _label, const size_t& _count)
{
	size_t quotient(_kmer / HTSIZE);
	htCell<HKMERr, ELMTr> e(quotient, _label, _count);
	return insert(_kmer, e);
}

//// used for loading without counts
	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::insert(const uint64_t& _kmer, const ILBL& _label)
{
	size_t q( _kmer / HTSIZE);
	htCell<HKMERr, ELMTr> e(q, _label, 1);
	size_t xline = _kmer % HTSIZE;
	m_table[xline].push_back(e);
	m_load++;
	return true;
}

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr,  ELMTr>::insert(const uint64_t& _kmer, const htCell<HKMERr, ELMTr>& _cell)
{
	size_t xline  = _kmer % HTSIZE;
	size_t q = _kmer  / HTSIZE;
	if (m_table[xline].empty())
	{
		//// min
		m_table[xline].push_back(_cell);
		//// max
		m_table[xline].push_back(_cell);
		//// cell iteself
		m_table[xline].push_back(_cell); 
	}
	else
	{
		//// insert
		m_table[xline].push_back(_cell);
		//// update min
		if (q < m_table[xline][0].CKey)
		{
			m_table[xline][0] = q;	
		}
		//// update max
		else if (q > m_table[xline][1].CKey)
		{
			m_table[xline][1] = q;
		}
	}
	m_load++;
	return true;
}

template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::find(const uint64_t& _kmer, size_t& _xElement, size_t& _yElement, ILBL& _label, IOCCR& _mult, ICount& _count) const
{
	const size_t remainder = _kmer % HTSIZE;
	if (m_table[remainder].empty())
	{
		return false;
	}
	HKMERr quotient = (_kmer / HTSIZE);
	if (quotient < m_table[remainder][0].CKey || quotient > m_table[remainder][1].CKey)
	{
		return false;
	}
	bool isContained = false;

	size_t t = 0;
	size_t t_index = 0;	

	size_t init_t =  rand() % (m_table[remainder].size() - 2);
	bool  dir_t = rand() % 2;

	for(t = init_t + 2; !isContained && t < m_table[remainder].size() && t >= 2; t = dir_t?t+1:t-1) 
	{
		isContained = m_table[remainder][t].CKey == quotient ;
		t_index = t;
	}
	if (isContained)
	{
		_xElement = remainder;
		_yElement = t_index;
		_label	  = m_table[remainder][t_index].CElement.GetLabel();
		_mult	  = m_table[remainder][t_index].CElement.GetMultiplicity();
		_count	  = m_table[remainder][t_index].CElement.GetCount();
		return true;
	}
	init_t = (init_t > 0 && init_t < m_table[remainder].size() - 2)? (dir_t ? init_t -1: init_t+1): init_t;

	for(t = init_t + 2;  !isContained && t < m_table[remainder].size() && t >= 2 ; t = dir_t?t-1:t+1)
	{
		isContained = m_table[remainder][t].CKey == quotient ;
		t_index = t;
	}
	if (isContained)
	{
		_xElement = remainder;
		_yElement = t_index;
		_label    = m_table[remainder][t_index].CElement.GetLabel();
		_mult     = m_table[remainder][t_index].CElement.GetMultiplicity();
		_count    = m_table[remainder][t_index].CElement.GetCount();
	}
	return isContained;
}

template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::find(const IKMER& _ikmer,  const size_t& _remainderI, const size_t& _quotientI, ILBL& _label) const
{
	const size_t remainder = _ikmer.skmer[_remainderI];

	if (m_table[remainder].empty())
	{
		return false;
	}
	size_t _endI = m_table[remainder].size() - 1;
	if (m_table[remainder][0].CKey > _ikmer.skmer[_quotientI] || m_table[remainder][_endI].CKey < _ikmer.skmer[_quotientI])
	{
		return false;
	}
	size_t midPoint;
	size_t _startI = 0;
	while (_endI > _startI)
	{
		midPoint = _startI + (_endI - _startI)/2;
		if (_ikmer.skmer[_quotientI] <= m_table[remainder][midPoint].CKey)
		{
			_endI = midPoint;
			continue;
		}
		_startI = midPoint+1;
	}
	if (m_table[remainder][_startI].CKey == _ikmer.skmer[_quotientI])
	{
		_label = m_table[remainder][_startI].CElement.Label;
		return true;
	}
	return false;
}

/* original function of the one below
	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::find(const uint64_t& _ikmer, ILBL& _label) 
{
	size_t quotient = _ikmer / HTSIZE;
	size_t remainder = _ikmer - quotient * HTSIZE;

	if (m_table[remainder].empty())
	{
		size_t _ikmerR = _ikmer;
		// The following 6 lines come from Jellyfish source code
		_ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
		_ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
		_ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
		_ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
		_ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
		_ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));

		quotient = _ikmerR / HTSIZE;
		remainder = _ikmerR - quotient * HTSIZE;

		if (m_table[remainder].empty())
		{       return false;   }
		uint8_t _endI = m_table[remainder].size() - 1;

		if (m_table[remainder][0].CKey > quotient || m_table[remainder][_endI].CKey < quotient)
		{       return false;}
		htCell<HKMERr, ELMTr>* ptr = &m_table[remainder].front();
		while (ptr->CKey <= quotient)
		{
			if (ptr->CKey == quotient)
			{
				_label = (ptr->CElement).Label;
				return true;
			}
			ptr++;
		}
		return false;
	}
	uint8_t _endI = m_table[remainder].size() - 1;
	if (m_table[remainder][0].CKey > quotient || m_table[remainder][_endI].CKey < quotient)
	{
		size_t _ikmerR = _ikmer;
		// The following 6 lines come from Jellyfish source code
		_ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
		_ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
		_ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
		_ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
		_ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
		_ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));

		quotient = _ikmerR / HTSIZE;
		remainder = _ikmerR - quotient * HTSIZE;

		if (m_table[remainder].empty())
		{       return false;   }
		_endI = m_table[remainder].size() - 1;
		if (m_table[remainder][0].CKey > quotient || m_table[remainder][_endI].CKey < quotient)
		{       return false;	}

		htCell<HKMERr, ELMTr>* ptr = &m_table[remainder].front();
		while (ptr->CKey <= quotient)
		{
			if (ptr->CKey == quotient)
			{
				_label = (ptr->CElement).Label;
				return true;
			}
			ptr++;
		}
		return false;
	}

	htCell<HKMERr, ELMTr>* ptr = &m_table[remainder].front();
	while (ptr->CKey <= quotient)
	{
		if (ptr->CKey == quotient)
		{
			_label = (ptr->CElement).Label;
			return true;
		}
		ptr++;
	}

	size_t _ikmerR = _ikmer;
	// The following 6 lines come from Jellyfish source code
	_ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
	_ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
	_ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
	_ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
	_ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
	_ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));

	quotient = _ikmerR / HTSIZE;
	remainder = _ikmerR - quotient * HTSIZE;
	//////////////////////////////////////////////////////////////////////////////////////////////////
	if (m_table[remainder].empty())
	{       return false;   }
	_endI = m_table[remainder].size() - 1; 
	if (m_table[remainder][0].CKey > quotient || m_table[remainder][_endI].CKey < quotient)
	{       return false;	}

	ptr = &m_table[remainder].front();
	while (ptr->CKey <= quotient)
	{
		if (ptr->CKey == quotient)
		{
			_label = (ptr->CElement).Label;
			return true;
		}
		ptr++;
	}
	return false;
	/////////////////////////////////////////////////////////////////////////////////////////////////////
}
*/

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::find(const uint64_t& _ikmer, ILBL& _label) 
{
	//// getting reverse kmer
	size_t _ikmerR = _ikmer;
	// The following 6 lines come from Jellyfish source code
	_ikmerR = ((_ikmerR >> 2)  & 0x3333333333333333UL) | ((_ikmerR & 0x3333333333333333UL) << 2);
	_ikmerR = ((_ikmerR >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_ikmerR & 0x0F0F0F0F0F0F0F0FUL) << 4);
	_ikmerR = ((_ikmerR >> 8)  & 0x00FF00FF00FF00FFUL) | ((_ikmerR & 0x00FF00FF00FF00FFUL) << 8);
	_ikmerR = ((_ikmerR >> 16) & 0x0000FFFF0000FFFFUL) | ((_ikmerR & 0x0000FFFF0000FFFFUL) << 16);
	_ikmerR = ( _ikmerR >> 32                        ) | (_ikmerR                        << 32);
	_ikmerR = (((uint64_t)-1) - _ikmerR) >> (64 - (m_k << 1));
	
	//// getting canonical kmer
	size_t _ikmerC = _ikmer < _ikmerR ? _ikmer : _ikmerR;

	size_t quotient = _ikmerC / HTSIZE;
	size_t remainder = _ikmerC - quotient * HTSIZE;

	//// check remainder in HT
	if (m_table[remainder].empty()) return false;
	
	//// check quotient range in HT
	uint8_t _endI = m_table[remainder].size() - 1;
	if (m_table[remainder][0].CKey > quotient || m_table[remainder][_endI].CKey < quotient)	return false;

	//// check quotients
	htCell<HKMERr, ELMTr>* ptr = &m_table[remainder].front();
	while (ptr->CKey <= quotient)
	{
		if (ptr->CKey == quotient)
		{
			_label = (ptr->CElement).Label;
			return true;
		}
		ptr++;
	}
	return false;
}

	template <typename HKMERr, typename ELMTr>
void hTable<HKMERr, ELMTr>::updateElement(const size_t& _xElement, const size_t& _yElement, const size_t& _count, const bool& _lbl, const bool& _isSameLbl)
{
	if (!_isSameLbl)
	{
		m_table[_xElement][_yElement].CElement.IncreaseMultiplicity();
	}
	m_table[_xElement][_yElement].CElement.AddToCount(_count);
	if (_lbl)
	{
		m_table[_xElement][_yElement].CElement.IncreaseMultiplicity();
	}
}

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::resetIterator()
{
	m_it_x = 0;
	m_it_y = 1;
	return true;
}

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::nextIterator()
{
	if (!m_table[m_it_x].empty() && m_table[m_it_x].size() > m_it_y + 1)
	{
		m_it_y++;
	}
	else
	{
		m_it_y = 2;
		m_it_x++;
	}
	return (!m_table[m_it_x].empty() && m_table[m_it_x].size() > m_it_y);
}

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::next()
{
	bool loop = true;
	while (loop && m_it_x < m_table.size())
	{
		loop = !nextIterator();
	}
	if (m_it_x >= m_table.size())
	{	return false;}
	return true;
}

	template <typename HKMERr, typename ELMTr>
void  hTable<HKMERr, ELMTr>::updateLabelAtIterator(const ILBL& _label)
{
	m_table[m_it_x][m_it_y].CElement.Label = _label;
}

	template <typename HKMERr, typename ELMTr>
void hTable<HKMERr, ELMTr>::markElementAtIterator()
{
	m_table[m_it_x][m_it_y].CElement.Mark();	
}

template <typename HKMERr, typename ELMTr>
void hTable<HKMERr, ELMTr>::elementIterator(ELMTr& _cElement) const
{
	_cElement = m_table[m_it_x][m_it_y].CElement;
}

template <typename HKMERr, typename ELMTr>
void hTable<HKMERr, ELMTr>::elementIterator(uint64_t& _kmer, ELMTr& _cElement) const
{
	_kmer = m_it_x + m_table[m_it_x][m_it_y].CKey*m_table.size();
	_cElement = m_table[m_it_x][m_it_y].CElement;
}

	template <typename HKMERr, typename ELMTr>
uint64_t hTable<HKMERr, ELMTr>::write(const char* _fileht, const size_t& _iteratorPos, const bool& _clearAfter)
{
	char * file_lbl = (char*) calloc(strlen(_fileht)+4,sizeof(char));
	char * file_key = (char*) calloc(strlen(_fileht)+4,sizeof(char));
	char * file_sze = (char*) calloc(strlen(_fileht)+4,sizeof(char));
	sprintf(file_lbl, "%s.lb", _fileht);
	sprintf(file_key, "%s.ky", _fileht);
	sprintf(file_sze, "%s.sz", _fileht);

	FILE * fd_l = fopen(file_lbl,"w+");
	FILE * fd_k = fopen(file_key,"w+");
	FILE * fd_s = fopen(file_sze,"w+");
	uint64_t nbElement = 0;
	uint8_t size = 0;
	for(ITYPE t = 0; t < HTSIZE; t++)
	{
		size = 0;
		if (!m_table[t].empty())
		{
			size_t  l_size = 0;
			//// count marked elements in bucket
			for(size_t u = _iteratorPos;  u < m_table[t].size() ;u++)
			{
				l_size += m_table[t][u].CElement.Marked() ? 1: 0;
			}
			if ( l_size >=  256) 
			{
				cerr << "This table can not be stored on disk: Some bucket list size exceeds 255." << endl;
				cerr << "Please relaunch all computations by applying the following modifications: " << endl;
				cerr << "- choose a smaller k-mers length, and/or" << endl;
				cerr << "- increase the size of the hash-table as " << m_table.size() << " is way too small!" << endl;
				cerr << "The program must exit now." << endl;
				exit(-1);
			}
			nbElement += l_size;
			size = l_size;
			//// write bucket size to file
			fwrite(&size, 1,1, fd_s);
			//// write keys and labels to files
			for(size_t u = _iteratorPos; u < m_table[t].size() ;u++)
			{
				if (m_table[t][u].CElement.Marked())
				{
					fwrite(&m_table[t][u].CKey,sizeof(HKMERr),1, fd_k);
					fwrite(&m_table[t][u].CElement.Label, sizeof(ILBL),1, fd_l);
				}
			}
			if (_clearAfter)
			{	
				m_table[t].clear();
			}
		}
		else
		{
			//// write bucket size 0 to file
			fwrite(&size, 1,1, fd_s);
		}
	}
	fclose(fd_l);
	fclose(fd_k);
	fclose(fd_s);
	free(file_lbl); 
	file_lbl=NULL;
	free(file_key); 
	file_key=NULL;
	free(file_sze); 
	file_sze=NULL;
	if (_clearAfter)
	{       
		m_table.clear();
	}
	return nbElement;
}

	template <typename HKMERr, typename ELMTr>
bool hTable<HKMERr, ELMTr>::read(const char * _filename, size_t& _fileSize,  const size_t& _nbCPU, const ITYPE& _modCollision, const bool& _isfastLoadingRequested)
{
	char * file_lbl = (char*) calloc(strlen(_filename)+4,sizeof(char));
	char * file_key = (char*) calloc(strlen(_filename)+4,sizeof(char));
	char * file_sze = (char*) calloc(strlen(_filename)+4,sizeof(char));

	sprintf(file_lbl, "%s.lb", _filename);
	sprintf(file_key, "%s.ky", _filename);
	sprintf(file_sze, "%s.sz", _filename);

	if (_isfastLoadingRequested)
	{
#ifdef _OPENMP
		omp_set_num_threads(_nbCPU);
#endif

		// Opening File with sizes
		_fileSize = HTSIZE;
		int fd_s = open(file_sze, O_RDONLY);
		if (fd_s == -1)
		{
			cerr << "Failed to open " << file_sze << endl;
			return false;
		}
		uint8_t *map;
		map = (uint8_t*) mmap(0, _fileSize, PROT_READ, MAP_SHARED, fd_s, 0);
		if (map == MAP_FAILED)
		{
			close(fd_s);
			cerr << "Failed to mmapping the file!" << endl;
			exit(-1);
		}
		if (_fileSize < 1)
		{
			close(fd_s);
			cerr << "Failed to load database: the database ["<< file_sze<< "] contains no data." << endl;
			exit(-1);
		}
		// Initialization
		ITYPE loadf = 0;
		bool allCollision = _modCollision <= 1;
		
		uint64_t 		nbElement = 0;
		vector<uint8_t>		choice(HTSIZE,0);
		vector<uint64_t> 	it_Key(_nbCPU,0);

		ITYPE i = 0;
		vector<ITYPE> Pos(_nbCPU+1, 0);
		for(i = 1; i < _nbCPU ; i++)
		{       Pos[i] = (HTSIZE/_nbCPU) * i;   }
		Pos[_nbCPU] = HTSIZE;

		/// PART 1: Setting bucket size
		i = 0;
		for(size_t t = 0; t < HTSIZE; t++)
		{
			if (t == Pos[i+1])
			{       it_Key[++i] = nbElement; }
			if (map[t] > 0)
			{
				loadf++;
				nbElement += map[t];
				choice[t] = (allCollision || (loadf % _modCollision)== 0) ? 2: 1;	
			}
		}
		i = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < _nbCPU; i++)
		{
			ITYPE min = Pos[i], max = Pos[i+1];
			for(ITYPE u = min; u < max; u++)
			{
				if (choice[u] == 2)
				{	m_table[u].resize(map[u]); } 
			}
		}
		/// PART 2: Populating key/label
		// Opening Files
		size_t _fileSizek = nbElement * sizeof(HKMERr);
		size_t _fileSizel = nbElement * sizeof(ILBL);
		// kmers keys
		int fd_k = open(file_key, O_RDONLY);
		if (fd_k == -1)
		{
			cerr << "Failed to open " << file_key << endl;
			return false;
		}
		HKMERr *key;
		key = (HKMERr*) mmap(0, _fileSizek, PROT_READ, MAP_SHARED, fd_k, 0);
		if (key == MAP_FAILED)
		{
			close(fd_k);
			cerr << "Failed to mmapping the file!" << endl;
			exit(-1);
		}
		if ( _fileSizek < 1)
		{
			close(fd_k);
			cerr << "Failed to load database: ["<< file_key<< "] contains no data." << endl;
			exit(-1);
		}
		/// targets Labels
		int fd_l = open(file_lbl, O_RDONLY);
		if (fd_l == -1)
		{
			cerr << "Failed to open " << file_lbl << endl;
			return false;
		}
		ILBL *lbl;
		lbl = (ILBL*) mmap(0, _fileSizel, PROT_READ, MAP_SHARED, fd_l, 0);
		if (lbl == MAP_FAILED)
		{
			close(fd_l);
			cerr << "Failed to mmapping the file!" << endl;
			exit(-1);
		}
		if ( _fileSizel < 1)
		{
			close(fd_l);
			cerr << "Failed to load database: ["<< file_lbl<< "] contains no data." << endl;
			exit(-1);
		}
		/// Loading in memory key-label for each k-mer
		i = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i = 0; i < _nbCPU ; i++)
		{
			ITYPE min = Pos[i], max = Pos[i+1];
			uint64_t u =0, v = 0, it_e = it_Key[i];

			for(ITYPE t = min; t < max; t++)
			{
				if (choice[t] > 0)
				{
					if (choice[t] == 2)
					{
						u = 0;
						for(v = it_e; v < it_e + map[t]; v++)
						{
							m_table[t][u].CKey = key[v];
							m_table[t][u++].CElement.Label = lbl[v];
						}
					}
					it_e += map[t];
				}
			}
		}
		/// End
		if (munmap(map, _fileSize) == -1)
		{
			perror("Error un-mmapping the database file (size)");
		}
		close(fd_s);
		if (munmap(key, _fileSizek) == -1)
		{
			perror("Error un-mmapping the database file (key)");
		}
		if (munmap(lbl, _fileSizel) == -1)
		{
			perror("Error un-mmapping the database file (label)");
		}
		close(fd_l);
		close(fd_k);

		_fileSize = HTSIZE + _fileSizek + _fileSizel;

		free(file_lbl); 
		file_lbl=NULL;
		free(file_key); 
		file_key=NULL;
		free(file_sze); 
		file_sze=NULL;

		return true;
	}
	FILE * fd_l = fopen(file_lbl,"r");
	FILE * fd_k = fopen(file_key,"r");
	FILE * fd_s = fopen(file_sze,"r");

#define LEN 100000
	
	if (fd_l == NULL)
	{	cerr << "Failed to open " << file_lbl << endl; return false;	}
	if (fd_k == NULL)
	{       cerr << "Failed to open " << file_key << endl; return false;   }
	if (fd_s == NULL)
	{       cerr << "Failed to open " << file_sze << endl; return false;   }


	ITYPE t = 0, i	= 0;
	ITYPE loadf 	= 0;
	bool allCollision = _modCollision <= 1;

	uint8_t 	c[LEN];
	ILBL 		lbl[LEN];
	HKMERr  	key[LEN];
	size_t len 	= fread(c,   1, LEN, fd_s);
	size_t len_l 	= fread(lbl, sizeof(ILBL), LEN, fd_l);
	size_t len_k 	= fread(key, sizeof(HKMERr), LEN, fd_k);

	size_t nbElement = 0, u = 0, v = 0, v_c = 0;

	_fileSize = len + len_l*sizeof(ILBL) + len_k*sizeof(HKMERr);
	//htCell<HKMERr, ELMTr> defCell;
	bool readDone = false;
	while (true)
	{
		while (i < len)
		{
			if (c[i] > 0)
			{
				loadf++;
				nbElement += c[i];
				readDone = false;
				if (allCollision || (loadf % _modCollision)== 0)
				{	
					//m_table[t].reserve(c[i]);	
					m_table[t].resize(c[i]);
					u = 0;
					v_c = v;
					while (u < c[i])
					{
						if (v_c < len_l)
						{
							//m_table[t].push_back(defCell);
							m_table[t][u].CKey = key[v_c];
							m_table[t][u++].CElement.Label = lbl[v_c++];
							continue;
						}
						len_l = fread(lbl, sizeof(ILBL), LEN, fd_l);
						len_k = fread(key, sizeof(HKMERr), LEN, fd_k);
						readDone = true;
						_fileSize += len_l*sizeof(ILBL) + len_k*sizeof(HKMERr);
						v_c = 0;
						if (len_l == 0)
						{       break;  }
					}
				}
				if (readDone)
				{	v = v_c;}
				else
				{
					if (v+c[i] >= LEN)
					{ 
						len_l = fread(lbl, sizeof(ILBL), LEN, fd_l);
						len_k = fread(key, sizeof(HKMERr), LEN, fd_k);
						_fileSize += len_l*sizeof(ILBL) + len_k*sizeof(HKMERr);
						if (len_l == 0)
						{       break;  }
						v = v+c[i] - LEN;
					}
					else
					{	v += c[i];}
				}  
			}
			i++;
			t++;
		}
		len = fread(c, 1, LEN, fd_s);
		_fileSize += len;
		i = 0;
		if (len == 0)
		{	break;	}
	}
	fclose(fd_l);
	fclose(fd_k);
	fclose(fd_s);

	free(file_lbl); 
	file_lbl=NULL;
	free(file_key); 
	file_key=NULL;
	free(file_sze); 
	file_sze=NULL;

	return true;	
}

