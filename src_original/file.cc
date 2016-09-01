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
#include <string>
#include <vector>
#include <string.h>

using namespace std;

#include "./file.hh"

void getElementsFromLine(char*& line, const size_t& len, const int _maxElement, std::vector< std::string >& _elements)
{
	size_t t = 0; 
	size_t cpt = 0;
	_elements.resize(0);
	while (t < len && cpt < _maxElement)
	{
		while ( t < len  && (line[t] == ' ' || line[t] == '\t' || line[t] == '\n' ))
		{       t++;
		}
		string v = "";
		while ( t < len && line[t] != ' '  && line[t] != '\t' && line[t] != '\n' )
		{
			v.push_back(line[t]);
			t++;
		}
		if (v != "")
		{
			_elements.push_back(v);
			cpt++;
		}
	}
	return;
}

void getElementsFromLine(const std::string& line, const size_t& _maxElement, std::vector< std::string >& _elements)
{
	size_t t = 0, len = line.size();
	size_t cpt = 0;
	_elements.resize(0);
	while (t < len && cpt < _maxElement)
	{
		while (t <  len && (line[t] == ' ' || line[t] == ',' || line[t] == '\n' || line[t] == '\t' ))
		{       t++;
		}
		string v = "";
		while (t <  len && line[t] != ' '  && line[t] != ',' && line[t] != '\n' && line[t] != '\t' )
		{
			v.push_back(line[t]);
			t++;
		}
		if (v != "")
		{	
			_elements.push_back(v);
			cpt++;
		}
	}
	return;
}

void getElementsFromLine(const std::string& line, const vector<char>& _seps, std::vector< std::string >& _elements)
{
	size_t t = 0, len = line.size();
	size_t cpt = 0;
	_elements.resize(0);
	while (t < len)
	{
		bool checkSep = true;
		while (t < len && checkSep)
		{
			checkSep = false;
			for(size_t i = 0; i < _seps.size() && !checkSep;i++)
			{ checkSep =  line[t] == _seps[i];}
			t += checkSep ? 1:0;
		}
		string v = "";
		checkSep = true;
		while (checkSep && t < len)
		{
			for(size_t i = 0 ; i < _seps.size() && checkSep; i++)
			{
				checkSep = checkSep && line[t] != _seps[i];
			}
			if (checkSep)
			{
				v.push_back(line[t]);
				t++;
			}
		}
		if (v != "")
		{_elements.push_back(v);}
		//cpt++;
	}
	return;
}

bool getLineFromFile(FILE*& _fileStream, string& _line)
{
	char *line = NULL;
	size_t len = 0;
	if (getline(&line, &len, _fileStream) != -1)
	{
		string l(line);
		if (l[l.size() -1 ] == '\n' )
		{l.erase(l.size()-1,1);}
		_line = l;
		free(line);
		line = NULL;
		return true;
	}
	else
	{
		//fclose(_fileStream);                       
		_line = "";
		free(line);
		line = NULL;
		return false;
	}
}

bool getFirstElementInLineFromFile(FILE*& _fileStream, string& _line)
{
	char *line = NULL;
	size_t len = 0;
	if (getline(&line, &len, _fileStream) != -1)
	{
		vector<string> ele;
		getElementsFromLine(line, len, 1, ele);
		_line = ele[0];
		free(line);
		line = NULL;
		return true;
	}
	else
	{
		_line = "";
		free(line);
		line = NULL;
		return false;
	}
}

bool getFirstAndSecondElementInLine(FILE*& _fileStream, uint64_t& _kIndex, ITYPE& _index)
{
	char *line = NULL;
	size_t len = 0;
	if (getline(&line, &len, _fileStream) != -1)
	{
		// Take the first element and put it into _kIndex: type IKMER
		// Take the second element and put it into _index: type ITYPE
		vector<string> ele;
		getElementsFromLine(line, len, 2, ele);

		_kIndex = atoll(ele[0].c_str());
		_index = atol(ele[1].c_str());
		free(line);
		line = NULL;
		return true;
	}
	else
	{
		free(line);
		line = NULL;
		return false;
	}
}

bool getFirstAndSecondElementInLine(FILE*& _fileStream, std::string& _line, ITYPE& _freq)
{
	char *line = NULL;
	size_t len = 0;
	if (getline(&line, &len, _fileStream) != -1)
	{
		// Take first element and put it into _line
		// Take second element and put it into _freq
		vector<string> ele;
		getElementsFromLine(line, len, 2, ele);
		_line = ele[0];
		_freq = atoi(ele[1].c_str());
		free(line);
		line = NULL;
		return true;
	}
	else
	{
		free(line);
		line = NULL;
		return false;
	}
}

