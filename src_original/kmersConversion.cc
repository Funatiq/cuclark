/*
 * CLARK, CLAssifier based on Reduced K-mers.
 */

/*
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

   Copyright 2013-2015, Rachid Ounit <rouni001@cs.ucr.edu>
 */

/*
 * @author: Rachid Ounit, Ph.D Candidate.
 * @project: CLARK, Metagenomics and Genomics Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
using namespace std;

#include "kmersConversion.hh"
#include "dataType.hh"

uint64_t getReverse(uint64_t _km_r, uint8_t k)
{
	// Code from Jellyfish
        _km_r = ((_km_r >> 2)  & 0x3333333333333333UL) | ((_km_r & 0x3333333333333333UL) << 2);
        _km_r = ((_km_r >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_r & 0x0F0F0F0F0F0F0F0FUL) << 4);
        _km_r = ((_km_r >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_r & 0x00FF00FF00FF00FFUL) << 8);
        _km_r = ((_km_r >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_r & 0x0000FFFF0000FFFFUL) << 16);
        _km_r = ( _km_r >> 32                        ) | ( _km_r                        << 32);
        _km_r = (((uint64_t)-1) - _km_r) >> (64 - (k << 1));
        return _km_r;
}

void getKmers(const std::string& c, uint64_t& _km_f, uint8_t k)
{
	_km_f = 0;
        uint8_t cpt = 0;
        while (cpt < k)
        {
                _km_f <<=2;
                switch (c[cpt++])
                {
                        case 'A': case 'a':_km_f ^= 0; break;
                        case 'C': case 'c':_km_f ^= 1; break;
                        case 'G': case 'g':_km_f ^= 2; break;
                        case 'T': case 't':_km_f ^= 3; break;
                        default:
                        cerr << "Failed to compute k-mer value of " << c << endl;
                        std::exit(1);
                        break;
                }
        }
}

void getReverseComplement(const uint64_t& _ikmer, const size_t& _kmerSize, uint64_t& _rev_kmerIndex)
{
	uint64_t tmp = _ikmer;
        _rev_kmerIndex = getReverse(tmp, (uint8_t)_kmerSize);
}

void getReverseComplement(const std::string& _kmer, uint64_t& _rev_kmerIndex)
{
	uint64_t tmp = 0;
        getKmers(_kmer, tmp,  (uint8_t)_kmer.size());
        _rev_kmerIndex = getReverse(tmp, (uint8_t)_kmer.size());
}

void vectorToIndex(const std::string& _kmer, uint64_t& _index)
{
	getKmers(_kmer, _index, _kmer.size());
}

void IndexTovector(const uint64_t& _index, const size_t& _kmerSize, std::string& _kmer)
{
        uint64_t num(_index);
        size_t max = 1;
        for(size_t t = 0 ; t < _kmerSize - 1; t++)
        {
                max *= 4;
        }
        vector<char> subKmer;
        while (max > 0)
        {
                if (num >= max)
                {
                        if (num/max == 1)
                        {
                                subKmer.push_back('C');
                                num = num - max;
                        }
                        if (num/max == 2)
                        {
                                subKmer.push_back('G');
                                num = num - 2*max;
                        }
                        if (num/max == 3)
                        {
                                subKmer.push_back('T');
                                num = num - 3*max;
                        }
                }
                else
                {
                        subKmer.push_back('A');
                }
                max = max / 4;
	}

        //Copy and reversing char order
        string k;
	for(size_t t = 0 ; t < subKmer.size() ; t++)
        {
                k.push_back(subKmer[subKmer.size() - 1 - t]);
      	} 
        _kmer = k;
}

