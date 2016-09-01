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

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <stdint.h>
#include "./file.hh"
#include <map>
using namespace std;

struct seqData
{
	std::string 	Name;
	uint64_t	GI;
	seqData():Name(""),GI(0)
	{}
};

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		cerr << "Usage: "<< argv[0] << " <./file of filenames> <./gi_taxid_nucl.dmp>"<< endl;
		exit(-1);
	}
	FILE * giToTx = fopen(argv[2], "r");
	if (giToTx == NULL)
	{
		cerr << "Failed to open " << argv[2] << endl;
		exit(-1);
	}
	FILE * meta_f = fopen(argv[1], "r");
	if (meta_f == NULL)
	{
		cerr << "Failed to open " << argv[1] << endl;
		exit(1);
	}
	///////////////////////////////////
	string line, file;
	vector<string> ele;
        vector<char> sep;
        sep.push_back('|');
	vector<int> TaxIDs;
	map<uint64_t,uint32_t>	 giToidx;
	vector<seqData> seqs;	
	map<uint64_t,uint32_t>::iterator it;
	uint32_t idx = 0;
	uint64_t gi = 0;
	cerr << "Loading GIs for all files... " ;
	while (getLineFromFile(meta_f, file))
	{
		FILE * fd = fopen(file.c_str(), "r");
		if (fd == NULL)
		{
			cerr << "Failed to open sequence file: " <<  file << endl;
			cout << file << "\tUNKNOWN" << endl;
			continue;
		}
		if (getLineFromFile(fd, line))
		{
			ele.clear();
			getElementsFromLine(line, sep, ele);
			if (ele.size() < 1 || ele[0] != ">gi")
			{
				continue;
			}
			gi = atol(ele[1].c_str());
			it = giToidx.find(gi);
			if (it == giToidx.end())
			{	
				TaxIDs.push_back(-1);
				giToidx[gi] = idx++;
			}
			seqData s;
			s.Name = file;
			s.GI = gi;
			seqs.push_back(s);
		}
		fclose(fd);
	}
	fclose(meta_f);
	cerr << "done" << endl;
	
	string pair;
	int taxID;
        vector<char> sepg;
        sepg.push_back(' ');
        sepg.push_back('\t');
	uint32_t cpt = 0, cpt_u = 0;
        cerr << "Retrieving taxonomy ID for each file... " ;
        size_t taxidTofind = TaxIDs.size(), taxidFound = 0;
	while (getLineFromFile(giToTx, pair) && taxidFound < taxidTofind)
        {
                ele.clear();
                getElementsFromLine(pair, sepg, ele);
                gi = atol(ele[0].c_str());
		taxID = atoi(ele[1].c_str());
                it = giToidx.find(gi);
		if (it != giToidx.end())
		{
			taxidFound++;
			TaxIDs[it->second] = taxID;
		}
        }
        fclose(giToTx);
	for(size_t t = 0; t < seqs.size(); t++)
	{
		cout << seqs[t].Name << "\t" ;
		it = giToidx.find(seqs[t].GI);
		cout << seqs[t].GI << "\t" << TaxIDs[it->second] << endl;
		if (TaxIDs[it->second]  == -1)
		{	cpt_u++; }
		else
		{	cpt++;	 }
	}
	cerr << "done (" << cpt << " files were successfully mapped";
	if (cpt_u > 0)
		cerr <<  ", "<< cpt_u << " files were unidentified" ;
	cerr << ")." << endl;
	return 0;
}

