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
	std::string	Accss; 
	seqData():Name(""),Accss("")
	{}
};

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cerr << "Usage: "<< argv[0] << " <./file of filenames> <./nucl_accession2taxid> <./merged.dmp>"<< endl;
		exit(-1);
	}
	FILE * oldTx = fopen(argv[3], "r");
        if (oldTx == NULL)
        {
                cerr << "Failed to open " << argv[3] << endl;
                exit(-1);
        }
	FILE * accToTx = fopen(argv[2], "r");
	if (accToTx == NULL)
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
	vector<string> ele, eles;
        vector<char> sep, seps;
        sep.push_back('|');
	sep.push_back('.');
	sep.push_back('>');
	seps.push_back(' ');
	seps.push_back('\t');
	seps.push_back(':');
	vector<int> TaxIDs;
	map<std::string,uint32_t> accToidx; 
	vector<seqData> seqs;	
	map<std::string,uint32_t>::iterator it;
	uint32_t idx = 0, i_accss = 0;
	string acc = "";
	cerr << "Loading accession number of all files... " ;
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
			getElementsFromLine(line, seps, ele);
			if (line[0] != '>' || ele.size() < 1)
			{
				continue;
			}
			eles.clear();
			getElementsFromLine(ele[0], sep, eles);
			
			i_accss = eles.size()>1?eles.size()-2:0;
			acc = eles[i_accss];
			it = accToidx.find(acc);
			
			if (it == accToidx.end())
			{	
				TaxIDs.push_back(-1);
				accToidx[acc] = idx++;
			}
			seqData s;
			s.Name = file;
			s.Accss = acc;
			seqs.push_back(s);
		}
		fclose(fd);
	}
	fclose(meta_f);
	cerr << "done ("<< accToidx.size() << ")" << endl;

	string on_line;
	sep.push_back(' ');
	sep.push_back('\t');
	std::map<int, int> 		oldTonew;
	std::map<int, int>::iterator	it_on;

	std::cerr << "Loading merged Tax ID... " ;
	while (getLineFromFile(oldTx,on_line))
	{
		ele.clear();
		getElementsFromLine(on_line, sep, ele);
		it_on = oldTonew.find(atoi(ele[0].c_str()));
		if (it_on == oldTonew.end())
		{
			oldTonew[atoi(ele[0].c_str())] = atoi(ele[1].c_str());
		}
	}	
	fclose(oldTx);
	std::cerr << "done" << std::endl;

	string pair;
	int taxID, new_taxID;
        vector<char> sepg;
        sepg.push_back(' ');
        sepg.push_back('\t');
	uint32_t cpt = 0, cpt_u = 0;
        cerr << "Retrieving taxonomy ID for each file... " ;
        size_t taxidTofind = TaxIDs.size(), taxidFound = 0;
	while (getLineFromFile(accToTx, pair) && taxidFound < taxidTofind)
        {
                ele.clear();
                getElementsFromLine(pair, sepg, ele);
                acc = ele[0];
		taxID = atoi(ele[2].c_str());
                it = accToidx.find(acc);
		if (it != accToidx.end())
		{
			taxidFound++;
			new_taxID = taxID;
			it_on = oldTonew.find(taxID);
			if (it_on != oldTonew.end())
			{	new_taxID = it_on->second;	}
			TaxIDs[it->second] = new_taxID;
		}
        }
        fclose(accToTx);
	for(size_t t = 0; t < seqs.size(); t++)
	{
		cout << seqs[t].Name << "\t" ;
		it = accToidx.find(seqs[t].Accss);
		cout << seqs[t].Accss << "\t" << TaxIDs[it->second] << endl;
		if (TaxIDs[it->second]  == -1)
		{	cpt_u++; }
		else
		{	cpt++;	 }
	}
	cerr << "done (" << cpt << " files were successfully mapped";
	if (cpt_u > 0)
	{	cerr <<  ", and "<< cpt_u << " unidentified";	}
	cerr << ")." << endl;
	return 0;
}

