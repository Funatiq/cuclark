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

#ifndef KMERSCONVERSION_HH
#define KMERSCONVERSION_HH

#include<string>
#include<stdint.h>

//      FUNCTION
//      Name: getReverseComplement
//      Implementation notes: Convert a k-mer number into its reverse k-mer.
//
void getReverseComplement(const uint64_t& _ikmer, const size_t& _kmerSize, uint64_t& _rev_kmerIndex);

void getReverseComplement(const std::string& _kmer, uint64_t& _rev_kmerIndex);

//      FUNCTION
//      Name: IndexTovector
//      Implementation notes: Convert a integer into a substring (or a sequence of bases).
//
void IndexTovector(const uint64_t& _index, const size_t& _kmerSize, std::string& _kmer);

//      FUNCTION
//      Name: vectorToIndex
//      Implementation notes: Convert a substring (or sequence of bases) into a 64-bit integer.
//
void vectorToIndex(const std::string& _kmer, uint64_t& _index);

#endif //KMERSCONVERSION_HH
