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
#include <map>
#include <stdint.h>
#include "./file.hh"
using namespace std;

#define NBNODE 6
struct node
{
	uint32_t 	parent;
	uint8_t		rank;
	node():parent(0),rank(255) {}
};

void getSGFOCP(const vector<node>& _nodes, const uint32_t& _taxid, vector<node>& _line)
{
	_line.clear();
	_line.resize(NBNODE);
	size_t cpt = 0;
	size_t it = _taxid, tmp;
	while (true)
	{
		if (it == 1 || _nodes[it].parent == 1)
			break;
		if (_nodes[it].rank < NBNODE && _line[ _nodes[it].rank ].rank != 0)
		{
			_line[ _nodes[it].rank ].rank = 0;
			_line[ _nodes[it].rank ].parent = it;
			cpt++;
		}
		tmp = it;
		it = _nodes[tmp].parent;
	}
}

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		cerr << "Usage: " << argv[0] << " <./nodes.dmp> <./file_taxid>"<< endl;
		exit(-1);
	}
	FILE * fdn = fopen(argv[1], "r");
	if (fdn == NULL)
	{
		cerr << "Failed to open " << argv[1] << endl;
		exit(-1);
	}
	FILE * fdt = fopen(argv[2], "r");
        if (fdt == NULL)
        {
                cerr << "Failed to open " << argv[2] << endl;
                exit(-1);
        }
	std::map<std::string,uint8_t> nameTorank;
	std::map<std::string,uint8_t>::iterator it;

	nameTorank["species"] = 0;
	nameTorank["genus"] = 1;
	nameTorank["family"] = 2;
	nameTorank["order"] = 3;
	nameTorank["class"] = 4;
	nameTorank["phylum"] = 5;

#define MAXNB 20000000
	vector<node> nodes(MAXNB);
	string line;
	vector<string> ele;
	vector<char> sep;
	sep.push_back(' ');
	sep.push_back('|');
	sep.push_back('\t');
	cerr << "Loading nodes of taxonomy tree... " ;
	int id, idp;
	while (getLineFromFile(fdn, line))
	{
		ele.clear();
		getElementsFromLine(line, sep, ele);
		id = atoi(ele[0].c_str());
		idp = atoi(ele[1].c_str());
		nodes[id].parent = idp;
		it = nameTorank.find(ele[2].c_str());
		if (it != nameTorank.end())
		{	nodes[id].rank = it->second; }
	}
	fclose(fdn);
	cerr << "done." << endl;
	vector<node> lineage;
	cerr << "Retrieving lineage for each sequence... " ;
	while (getLineFromFile(fdt, line))
	{
		ele.clear();
		getElementsFromLine(line, sep, ele);
		id = atoi(ele[2].c_str());
		cout << ele[0] << "\t" << id;
		if (id > 0)
		{	
			getSGFOCP(nodes, id, lineage);	
			for(size_t t = 0; t < NBNODE; t++)
			{
				if (lineage[t].rank == 0)
				{
					cout << "\t" << lineage[t].parent;
				}
				else
				{
					cout << "\tUNKNOWN";
				}
			}
			cout << endl;
			continue;
		}
		for(size_t t = 0; t < NBNODE; t++)
		{	cout << "\tUNKNOWN" ;
		}
		cout << endl;
	}
	fclose(fdt);
	cerr << "done." << endl;
	return 0;
}

