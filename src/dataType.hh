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
 * Added new types.
 */

#ifndef DATATYPE_HH
#define DATATYPE_HH

#include <stdint.h>
#include "./parameters.hh"
#include <iostream>
typedef uint32_t        ITYPE;
typedef uint8_t        	IOCCR;
typedef uint16_t	ILBL;

////// new types
typedef uint16_t		RESULTS;
//~ typedef uint8_t		CONTAINER;
typedef uint16_t		CONTAINER;
//~ typedef uint32_t		CONTAINER;
//////

#define MTRGTS 		65535

struct IKMER
{
	uint64_t         skmer[SB];

	IKMER() 
	{}

	IKMER(const uint64_t& _kmer)
	{
		uint64_t quotient(_kmer);
		uint64_t k = 0;
		while (k < SB)
		{
			skmer[k] = quotient % HTSIZE;
			quotient =  (quotient) / HTSIZE;
			k++;
		}
	}

	IKMER(const std::string& _kmer)
	{
                uint64_t sum = 0, sum_r = 0;
                size_t val;
                size_t k = _kmer.size();
                size_t r_t = 0;
                char l;
                for (size_t t = 0; t < k ; t++)
                {
                        r_t = k - t - 1;
                        l = _kmer[r_t];
                        val = l == 'A' ? 0:(l == 'C' ? 1:(l == 'G' ? 2:3));
                        sum = sum*4 + val;
                        l = _kmer[t];
                        val = l == 'A' ? 3:(l == 'C' ? 2:(l == 'G' ? 1:0));
                        sum_r = sum_r*4 + val;
                }
                skmer[0] = sum % HTSIZE;
                skmer[1] = sum/ HTSIZE;
                skmer[2] = sum_r % HTSIZE;
                skmer[3] = sum_r/ HTSIZE;
	}

	uint64_t getIKMER() const
	{
		uint64_t value = 0;
		uint64_t base = 1;
		for(size_t t = 0; t < SB; t++)
		{
			value += skmer[t]*base;
			base *= HTSIZE;
		}
		return (value);
	}

	bool operator==(const uint64_t& a) const
	{
		return (this->getIKMER() == a);
	}

	IKMER& operator=(const std::string& _kmer)
	{
		uint64_t sum = 0, sum_r = 0;
		size_t val;
		size_t k = _kmer.size();
		size_t r_t = 0;
		char l;
		for (size_t t = 0; t < k ; t++)
        	{
			r_t = k - t - 1;
                	l = _kmer[r_t];
			val = l == 'A' ? 0:(l == 'C' ? 1:(l == 'G' ? 2:3));
                	sum = sum*4 + val;
			l = _kmer[t];
			val = l == 'A' ? 3:(l == 'C' ? 2:(l == 'G' ? 1:0));
			sum_r = sum_r*4 + val;
        	}
		skmer[0] = sum % HTSIZE;
		skmer[1] = sum/ HTSIZE;
		skmer[2] = sum_r % HTSIZE;
		skmer[3] = sum_r/ HTSIZE;

		return *this;
	}

	void SetReverse(const uint64_t& _ikmer)
	{
		skmer[2] = _ikmer % HTSIZE;
		skmer[3] = _ikmer / HTSIZE;
	}

	IKMER& operator=(const uint64_t& _q)
	{
		size_t quotient(_q);
		size_t k = 0;
		while (k < SB)
		{
			skmer[k] = quotient % HTSIZE;
			quotient =  (quotient) / HTSIZE;
			k++;
		}
		return *this;
	}

};

struct ICount
{
	uint16_t         digit_0;
	uint16_t         digit_1;

	ICount(): digit_0(0), digit_1(0)
	{}

	ICount(const size_t& _count)
	{
		if (_count > 4294967296)
		{
			digit_0 = 65535;
			digit_1 = 65535;
		}
		else
		{
			digit_0 = _count % 65536;
			digit_1 = ((_count)/65536) % 65536;
		}
	}

	ICount& operator =(const size_t& _count)
	{
		if (_count > 4294967296)
		{
			digit_0 = 65535;
			digit_1 = 65335;
		}
		else
		{
			digit_0 = _count % 65536;
			digit_1 = ((_count)/65536) % 65536;
		}
		return *this;
	}

	bool operator==(const ICount& _count) const
	{
		return (digit_0 == _count.digit_0 && digit_1 == _count.digit_1);
	}

	bool operator >(const ICount& _count) const
	{
		return (this->getCount() > _count.getCount());
	}

	bool operator >=(const ICount& _count) const
	{
		return (this->getCount() >= _count.getCount()); 
	}

	size_t getCount() const
	{
		return (digit_0 + 65536*digit_1);
	}

};

template <typename T> struct sVector
{
        T*       ptr;
        uint32_t idx;
        uint32_t cap;

        sVector(): idx(0), cap(0), ptr(NULL)
        {}
        ~sVector()
        {
                clear();
        }
        size_t  size()  const
        {
                return idx;
        }
        bool empty() const
        {
                return (idx == 0);
        }
        T& front() const
        {
                return (*ptr);
        }
        void resize(const size_t& _size)
        {
                ptr = (T*) malloc(_size * sizeof(T));
                idx = _size;
                cap = _size;
        }
        void push_back(const T& _cell)
        {
                if (idx < cap)
                {
                        ptr[idx++] = _cell;
                }
                else
                {
                        // Creation of a new array
                        cap = cap == 0 ? 1: 2*cap;
                        T* ptr_n = (T*) malloc(cap * sizeof(T));
                        // Copy & update
                        for(size_t t = 0; t < idx; t++)
                        {
                                ptr_n[t] = ptr[t];
                        }
                        ptr_n[idx++] = _cell;
                        free(ptr);
                        ptr = NULL;
                        // Re-addressing
                        ptr = ptr_n;
                }
        }
        T& operator[](const size_t& _pos) const
        {
                return ptr[_pos];
        }
        T* begin() const
        {
                return ptr;
        }
        T* end() const
        {
                return ptr+idx;
        }
        void clear()
        {
                free(ptr);
                ptr = NULL;
                idx = 0;
                cap = 0;
        }
};

struct Element
{
	ILBL		Label;
	IOCCR           Multiplicity;
	ICount          Count;
	Element(): Label(0), Multiplicity(1), Count(){}

	void Set(const ILBL& _label, const size_t& _count) 
	{	Label = _label; 
		Count = _count; }
	void SetLabel(const ILBL& _label)
	{	Label = _label;}
	ILBL GetLabel() const	
	{	return Label;	}
	IOCCR GetMultiplicity() const 
	{	return Multiplicity;	}
	size_t GetCount() const 
	{	return Count.getCount();	}
	void AddToCount(const size_t& _count) 
	{	Count = Count.getCount() + _count; 	}
	void IncreaseMultiplicity() 
	{ 	Multiplicity += Multiplicity < 255? 1: 0;	}
	void Mark()	
	{	Count = 0;	}	
	bool Marked() const 
	{	return Count == 0;	}

};

struct lElement
{
        uint16_t	Label;
        uint8_t       	Multiplicity;
        uint8_t		Count;
	lElement(): Label(0), Multiplicity(1), Count(0)
	{}
        void Set(const ILBL& _label, const size_t& _count)
        {
		Label = _label;       
		Count = (uint8_t) _count;
	}
        ILBL GetLabel() const
        {       return Label;   }
        IOCCR GetMultiplicity() const
        {       return Multiplicity;    }
        size_t GetCount() const
        {       return ((size_t) Count);        }
        void AddToCount(const size_t& _count)
        { 	Count += ((size_t) Count) + _count < 255 ? _count: 0;	}
        void IncreaseMultiplicity()
        {       Multiplicity += Multiplicity < 255? 1: 0;       }
        void Mark()
        {       Count = 0; 	}	
        bool Marked() const
        {       return (Count == 0); 	}
};

struct rElement
{
        ILBL            Label;
        rElement(): Label(0){}
        void Set(const ILBL& _label, const size_t& _count)
        {
                Label = _label;
        }
	void SetLabel(const ILBL& _label)
	{	Label = _label;	}
        ILBL GetLabel() const {return Label;}
        IOCCR GetMultiplicity() const {return 1;}
        size_t GetCount() const { return 1;}
        void AddToCount(const size_t& _count) {}
        void IncreaseMultiplicity() {}
	void Mark()	{;}
	bool Marked() const {return true;}

};

struct ObjectData
{
	bool		BumpFound;
	int		MinCount;
	int		MaxCount;

	ObjectData(): BumpFound(false), MinCount(1), MaxCount(2000000) {}
};

#endif
