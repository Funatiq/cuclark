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
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>

#include "file.hh"

using namespace std;

#include "analyser.hh"

analyser::~analyser()
{
	m_kmers.clear();
	m_frequency.clear();
}
bool analyser::getBumpInterval(int& _indexS, int& _indexE, const size_t& _div)
{
	int unchanged = 0;
	vector< vector<int> > _freqTable;
	int _minFreq = 1, _maxFreq = 1, indexS = 0;
	int _minVal = 999999999, _maxVal = 0;
	bool minfound = false;
	for(size_t i =0 ; i < m_frequency.size() ; i++)
	{
		if (i > 0 && labs(m_frequency[i] - m_frequency[i-1]) < 1 )
		{
			if (unchanged < 5)
			{
				unchanged++;
			}
		}
		else
		{       unchanged = 0;

		}

		if (unchanged < 1)
		{
			vector<int> point(2);
			point[0] = i;
			point[1] = m_frequency[i];
			_freqTable.push_back(point);
			if (!minfound && point[1] > 0)
			{
				indexS = i;
				minfound = true;
			}
		}
	}
	size_t Length = _freqTable.size();
	if (Length < 3)
	{
		_indexS = indexS;
		_indexE = _freqTable[Length - 1][0];
		return false;
	}
	_indexS = indexS; 
	_indexE = _freqTable[Length - 1][0];
	bool ret = false;

	if (Length <= 4)
	{
		return ret;
	}
	else
	{
		bool minDone = false;
		bool maxDone = false;
		for(size_t cpt = 1; cpt < Length && !maxDone ; cpt++)
		{ 
			vector<int> point = _freqTable[cpt];
			if (!minDone && _minVal >= point[1])
			{
				_minVal = point[1];
				_minFreq = point[0];
			}

			int step = (0.5 * _minFreq) >= 2 ? 0.5 * _minFreq : 2;
			minDone = point[0] - _minFreq >= step; 

			if (!maxDone && minDone && _maxVal < point[1])
			{
				_maxVal = point[1];
				_maxFreq = point[0];
			} 

			maxDone = minDone && (point[0] - _maxFreq >= (_maxFreq - _minFreq));

		}
		if (maxDone)
		{
			_indexS = _maxFreq - (_maxFreq - _minFreq)/_div;
			_indexE = _maxFreq + (_maxFreq - _minFreq)/_div;
			ret =  true;
		}
	}
	for(size_t t = 0; t < Length ; t++)
	{
		_freqTable[t].clear();
	}
	_freqTable.clear();

	return ret;

}

analyser::analyser(const char* _file)
{
	m_size = 0;

	string n_s;
	int n_i = 0;
	int n_max = 0;

	ifstream n_infile(_file);
	while (n_infile >> n_s >> n_i)
	{
		if (n_i > n_max)
		{
			n_max = n_i;
		}
	}
	m_size = n_max + 1;
	n_infile.close();

	m_kmers.resize(m_size);
	m_frequency.resize(m_size);

	ifstream infile(_file);
	int cpt = 0;
	string s;
	int i = 0, index = 0;

	while (infile >> s >> i)
	{		
		m_frequency[i]++;
		m_kmers[i].push_back(index);		
		index++;
	}
	infile.close();
}

