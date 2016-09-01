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

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <fstream>
using namespace std;
#include "./file.hh"

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cerr << "Usage: " << argv[0] << " <FilestoTaxIDs>, option: <Rank: 0,1,2,3,4,5>, 0 for species,...,5 for phylum. Default is species." << endl; 
		exit(1);
	}
	FILE * fd = fopen(argv[1], "r");
	if (fd == NULL)
	{
		cerr << "Failed to open " << argv[1] << endl;
		exit(1);
	}
	int r = 1;
	if (argc > 2)
	{
		r = atoi(argv[2]);
		if (r > 5)
		{
			cerr << "Failed to recognize the rank. Please type a number between 0 and 5, according to the following:" << endl;
			cerr << "0: species, 1: genus, 2: family, 3: order, 4:class, and 5: phylum." << endl;
			exit(1);
		}
	}
	string line;
	vector<string> ele;
	vector<char> sep;
	sep.push_back('\t');
	sep.push_back(',');
	sep.push_back(' ');
	size_t inc = 0;
	size_t nbFilesExcluded = 0;
	ofstream fout("files_excluded.txt", std::ios::binary);
	while (getLineFromFile(fd, line))
	{
		ele.clear();
		getElementsFromLine(line, sep, ele);
		if ( ele[1] != "-1" )
		{
			if (ele[2+r] != "UNKNOWN")
			{
			cout << ele[0] << "\t" << ele[2+r] << endl;
			}
		}
		else
		{
			nbFilesExcluded++;
			if (nbFilesExcluded == 1)
			{ 
				fout << "The following files have been excluded from the targets definition" << endl;
			}
			// Report this file
			fout << ele[0] << endl;
		}
	}
	fclose(fd);
	fout.close();
	return nbFilesExcluded;
}
