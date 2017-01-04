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

   Copyright 2013-2016, Rachid Ounit <rouni001@cs.ucr.edu>
 */

/*
 * @author: Rachid Ounit, Ph.D Candidate.
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#ifndef ANALYSER_HH
#define ANALYSER_HH

/*
 * Rachid Ounit
 * Date: 07/12/2013
 * 
 */

#include <vector>

class analyser
{
	private:
		size_t				m_size;
		std::vector< std::vector<int> > m_kmers;
		std::vector<int> 		m_frequency;

	public:
		analyser(const char* _file);
		~analyser();
		bool getBumpInterval(int& _indexS, int& _indexE, const size_t& _div = 2);
		
};


#endif
