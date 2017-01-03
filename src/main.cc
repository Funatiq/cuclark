/*
   CuCLARK, CLARK for CUDA-enabled GPUs.
   Copyright 2016-2017, Robin Kobus <rkobus@students.uni-mainz.de>
   
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
 * Deprecated some parameters and added new parameters.
 */

#include<iostream>
#include<cstdlib>
#include<stdio.h>
#include<cmath>
#include<fstream>

#include "./CuCLARK_hh.hh"
#include "./parameters.hh"
#define MAXK 32

using namespace std;

void printUsage()	
{
	cout << "\n";
	cout << "CuCLARK -- ''CLARK for CUDA-enabled GPUs'' \n\n";
	
	cout << "./CuCLARK -k <kmerSize> -t <minFreqTarget> -T <fileTargets> -D <directoryDB/> -O <fileObjects> -R <fileResults> -n <numberofthreads> ...\n\n\n";
	cout << "Definitions of parameters (cf. README file, for details):\n\n";
	
	cout << "-k <kmerSize>,       \t k-mer length:\tinteger, >= 2 and <= 32. The default value is 31.\n";
	cout << "-t <minFreqTarget>,  \t minimum of k-mer frequency in targets:\tinteger, >=0.\n";
	//~ cout << "-o <minFreqtObject>, \t minimum of k-mer frequency in objects:\tinteger, >=0.\n";
	cout << "-T <fileTargets>,    \t filename of the targets definition:\t text.\n";
	cout << "-D <directoryDB/>,   \t filename of the database directory (to load/save database files):\t text.\n";
	cout << "-O <fileObjects>,    \t filename of objects (or list of objects):\t text.\n";
	cout << "-P <file1> <file2>,  \t filenames of paired-end reads:\t texts.\n";
	cout << "-R <fileResults>,    \t filename to store results (or corresponding list of results file):\t text.\n";
	//~ cout << "-m <mode>,           \t mode of execution: 0 (full), 1 (default), 2 (express) or 3 (spectrum).\n";
	cout << "-n <numberofthreads>,\t number of threads:\tinteger >= 1.\n";
	cout << "-b <numberofbatches>,\t number of batches:\tinteger >= 1.\n";
	cout << "-d <numberofdevices>,\t number of CUDA devices to use:\tinteger >= 1.\n";
	cout << "--tsk,               \t to request a detailed creation of the database (target specific k-mers files).\n";
	//~ cout << "--ldm,               \t to request the loading of the database by memory mapped-file (in multithreaded mode, multiple parallel threads are requested).\n";
	cout << "-g <iteration>,      \t gap or number of non-overlapping k-mers to pass for the database creation (for CuCLARK-l only). The default value is 4.\n";
	cout << "-s <factor>,         \t sampling factor value (for CuCLARK only).\n";
	cout << "\n";
	cout << "--help,              \t to print help/options.\n";
	cout << "--version,           \t to print the version info.\n";
	cout << endl;
	return;
}

int main(int argc, char** argv)
{
	// check number of arguments
	if (argc == 2)
	{
		string val(argv[1]);
		if (val == "--help" || val == "--HELP" )
		{       printUsage();   return 0; }
		if (val == "--version" || val == "--VERSION")
		{
			cout << "Version: " << VERSION << " (Copyright 2016-2017 Robin Kobus, rkobus@students.uni-mainz.de)" << endl;
			cout << "Based on CLARK version 1.1.3 (UCR CS&E. Copyright 2013-2016 Rachid Ounit, rouni001@cs.ucr.edu) " << endl;	
			return 0; 
		}
	}
	if (argc < 6)
	{		
		cerr << "To run " << argv[0] << ", at least four  parameters are necessary:\n" ;
		cerr << "filename of the targets definition, directory of database, filename for objects, filename for results."<< endl;
		printUsage();
		return -1;
	}
	
	// set default parameters
	size_t k = 31, mode = 0, cpu = 1, iterKmers = 0;
	ITYPE minT = 0, minO = 0, sfactor = 1;
	bool cLightDB = false, ldm = false, tsk = false, kso= false, ext = false, isReduced = false;
	int i_targets = -1, i_objects = -1, i_objects2 = -1, i_folder=-1, i_results =-1;
	char * mergedFiles = NULL;
	
	size_t batches = 1, dbParts = 1, devices = 0;

	// parse arguments
	for(size_t i = 1; i < argc ; i++)
	{
		string val(argv[i]);
		if (val ==  "-k")
		{	if (i++ >= argc) {cerr << "Please specify the k-mer length!"<< endl; exit(1);	}
			k =  atoi(argv[i]); 
			if (k <= 1 || k > MAXK) { cerr <<"The k-mer length should be in [2,"<< MAXK << "]." << endl; exit(1);}
			continue;
		}
		if (val ==  "-t")
		{
			if (i++ >= argc) {cerr << "Please specify the minimum frequency (targets)!"<< endl; exit(1);    }
			minT =  atoi(argv[i]);
			if (minT >= 65536) { cerr <<"The min k-mer frequency should be in [0,65535]." << endl; exit(1);}
			continue;
		}
		//// deprecated parameters
		//~ if (val ==  "-o")
		//~ {
			//~ if (i++ >= argc) {cerr << "Please specify the minimum frequency (objects)!"<< endl; exit(1);    }
			//~ minO =  atoi(argv[i]);
			//~ if ( minO >= 65536) { cerr <<"The min k-mer frequency should be in [0,65535]."  << endl; exit(1);}
			//~ continue;
		//~ }
		//~ if (val ==  "-m")
		//~ {
			//~ if (i++ >= argc) {cerr << "Please specify the mode!"<< endl; exit(1);    }
			//~ mode =  atoi(argv[i]);
			//~ if (mode > 3) {	cerr <<"The mode of execution should be 0, 1, 2, or 3." << endl; exit(1);}
			//~ if (mode != 3) {	kso = false;	}
			//~ if (mode != 1) {	sfactor = 1;	}
			//~ continue;
		//~ }
		////
		if (val ==   "-n")
		{
			if (i++ >= argc) {cerr << "Please specify the number of threads!"<< endl; exit(1);    }
			cpu =  atoi(argv[i]);
			if (batches < cpu) batches = cpu;
			if (cpu< 1) { cerr <<"The number of threads should be higher than 0." << endl; exit(1);}
			continue;
		}
		//// deprecated parameter
		//~ if (val ==   "--ldm")
		//~ {
			//~ ldm = true; continue;
		//~ }
		////
		if (val ==  "--tsk")
		{
			tsk = true; continue;
		}
		if (val ==   "--extended")
        {
            ext = true; continue;
        }
		if (val ==  "-T")
		{
			if (i++ >= argc) {cerr << "Please specify the targets!"<< endl; exit(1);    }
			i_targets =  i;
			if (!validFile(argv[i])) { cerr <<"Failed to find/read the file of the targets definition: " << argv[i] << endl; exit(1);}
			continue;
		}
		if (val ==  "-O")
		{
			if (i++ >= argc) {cerr << "Please specify the objects!"<< endl; exit(1);    }
			i_objects =  i;
			if (!validFile(argv[i])) { cerr <<"Failed to find/read the filename of objects: " << argv[i] << endl; exit(1);}
			continue;
		}
		if (val ==  "-P")
        {
            if (i+2 >= argc) {cerr << "Please specify the paired-end reads!"<< endl; exit(1);    }
            i++;
			i_objects =  i;
			i_objects2 = i+1;
            if (!validFile(argv[i++])) { cerr <<"Failed to find/read " << argv[i-1] << endl; exit(1);}
			if (!validFile(argv[i]))   { cerr <<"Failed to find/read " << argv[i] << endl; exit(1);}
			//~ mergedFiles = (char *) calloc(strlen(argv[i-1])+25, sizeof(char));
			//~ sprintf(mergedFiles,"%s_ConcatenatedByCLARK.fa",argv[i-1]);
			//~ mergePairedFiles(argv[i-1], argv[i], mergedFiles);
 			continue;
 		}
		if (val ==   "-D")
		{
			if (i++ >= argc) {cerr << "Please specify the database directory!"<< endl; exit(1);    }
			i_folder =  i;
			if (!validFile(argv[i])) { cerr <<"Failed to find/read the directory:  " << argv[i] << endl; exit(1);}
			continue;
		}
		if (val ==   "-R")
		{
			if (i++ >= argc) {cerr << "Please specify where to store results!"<< endl; exit(1);    }
			i_results =  i;
			continue;
		}
		if (val == "-g")
		{
			if (i++ >= argc) {cerr << "Please specify a gap value!"<< endl; exit(1);    }
                  iterKmers =  atoi(argv[i]);
			if (iterKmers < 4) {cerr << "The gap value should be >= 4."<< endl; exit(1);    }
            continue;
        }
		if (val == "-s")
        {
            if (i++ >= argc) {cerr << "Please specify a sampling factor value!"<< endl; exit(1);    }
            sfactor =  atoi(argv[i]);
            if (sfactor < 2 || sfactor > SFACTORMAX) 
			{	cerr << "The sampling factor value should be in the interval [2,"<< SFACTORMAX <<"]."<< endl; exit(1);    }
            continue;
        }
        // new parameters      
        if (val == "-b")
        {
			if (i++ >= argc) {cerr << "Please specify the number of batches!"<< endl; exit(1);    }
			batches =  atoi(argv[i]);
			if (batches < cpu)
			{	cerr << "The number of batches should be higher than the number of threads."<< endl; exit(1);    }
			continue;
		}
 		if (val == "-d")
        {
			if (i++ >= argc) {cerr << "Please specify the number of devices to use!"<< endl; exit(1);    }
			devices =  atoi(argv[i]);
			if (devices < 1)
			{	cerr << "The number of devices should be higher than 0."<< endl; exit(1);    }
			continue;
		}
		//
		           
		cerr << "Failed to recognize option: " << val << endl;
		exit(1);
	}
	
	if (HTSIZE == LHTSIZE)
	{	// CuCLARK-l
		cLightDB = true; 
		if (iterKmers == 0)
		{ 	iterKmers = 4 ;}
	 	k = 27;
		if (cLightDB)
		{	sfactor = 1;}
	}
	else
	{
		iterKmers = 0;
	}
	if (k == 0) // should never occur because of the code above
	{
		cerr << "Please specify a k-mers length: -k <integer>" << endl;
		exit(1);
	}
	// check for sufficient parameters
	if ( i_targets < 0 || i_folder < 0 || i_objects < 0 ||  i_results < 0)
	{
		cerr << "Failed to run " << argv[0] << ": at least four  parameters are necessary" ;
		cerr << ": file of targets, directory of database, file of objects, file for results."<< endl;

		printUsage();
		exit(1);
	}
	const char * objects 	= i_objects > 0 ? argv[i_objects] : NULL;
	const char * objects2 	= i_objects2 > 0 ? argv[i_objects2] : NULL;
	const bool paired 	= i_objects2 > 0;
	
	string folder(argv[i_folder]);
	if (folder[folder.size()-1] != '/')
	{	folder.push_back('/');	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	size_t t_b = log(HTSIZE)/log(4.0);	// 15 for (Cu)CLARK, 12 for (Cu)CLARK-l
	size_t max16 = t_b + 8;				// 23 for (Cu)CLARK, 20 for (Cu)CLARK-l
	size_t max32 = t_b + 16;			// 31 for (Cu)CLARK, 28 for (Cu)CLARK-l
	
	// (Cu)CLARK-l should have k == 27 and therefore max16 == 20 <= k <= 28 == max32

	if (k <= max16)
	{
		// Use 2Bytes to store each discriminative k-mer
		CuCLARK<T16> classifier(k, argv[i_targets], folder.c_str(), minT, tsk, cLightDB, iterKmers, cpu, sfactor, ldm, batches, devices);
		if (paired)
			classifier.run(objects, objects2, argv[i_results], minO, ext);
		else
			classifier.run(objects, argv[i_results], minO, ext);
		//~ deleteFile(mergedFiles);
		exit(0);
	}
	if (k <= max32)
	{
		// Use 4Bytes to store each discriminative k-mer
		CuCLARK<T32> classifier(k, argv[i_targets], folder.c_str(), minT, tsk, cLightDB, iterKmers, cpu, sfactor, ldm, batches, devices);
		if (paired)
			classifier.run(objects, objects2, argv[i_results], minO, ext);
		else
			classifier.run(objects, argv[i_results], minO, ext);
		//~ deleteFile(mergedFiles);
		exit(0);
	}
	if (k <= MAXK)
	{
		// Use 8Bytes to store each discriminative k-mer
		CuCLARK<T64> classifier(k, argv[i_targets], folder.c_str(), minT, tsk, cLightDB, iterKmers, cpu, sfactor, ldm, batches, devices);
		if (paired)
			classifier.run(objects, objects2, argv[i_results], minO, ext);
		else
			classifier.run(objects, argv[i_results], minO, ext);
		//~ deleteFile(mergedFiles);
		exit(0);
	}
	std::cout <<"This version of CuCLARK does not support k-mer length strictly higher than " << MAXK << std::endl;
	exit(-1);

}

