# CuCLARK

ABOUT
-----
CuCLARK is a metagenomic classifier for CUDA-enabled GPUs, based on CLARK (http://clark.cs.ucr.edu/).  
For implementation details and speed comparison see the corresponding paper [Accelerating metagenomic read classification on CUDA-enabled GPUs](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1434-6). CuCLARK [v1.0](https://github.com/Funatiq/cuclark/releases/tag/v1.0) was used in the paper and has since been updated (see `CHANGELOG.md` for details).


The program comes in two variants: CuCLARK and CuCLARK-l.
CuCLARK is designed for workstations which can provide enough RAM to fit large databases 
(like all bacterial genomes in NCBI/RefSeq) and an arbitrary number of CUDA-enabled GPUs.
CuCLARK-l uses a much smaller database and can be run on computers with 4 GB of RAM and 
a CUDA-enabled GPU with at least 1 GB of memory.

To use CuCLARK as a metagenome classifier, it is recommended to use the provided scripts 
for that purpose, which are detailed in the section "CLASSIFICATION OF METAGENOMIC SAMPLES".


AUTHOR
-----
Robin Kobus, master student at Institute of Computer Science, JGU Mainz, Germany


ORIGINAL CLARK AUTHORS
-----
Rachid Ounit (1), Steve Wanamaker(2), Timothy J Close (2) and
Stefano Lonardi (1). 
1: Computer Science and Engineering department, University of California, 
Riverside CA 92521
2: Botany and Plant Sciences department, University of California, 
Riverside CA 92521

Project dates: 08/2013 to 03/2015.


LICENSE
-----
CuCLARK, CLARK for CUDA-enabled GPUs.
Copyright 2016, Robin Kobus <rkobus@students.uni-mainz.de>
   
based on CLARK version 1.1.3, CLAssifier based on Reduced K-mers.
Copyright 2013-2016, Rachid Ounit <rouni001@cs.ucr.edu>

CuCLARK is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CuCLARK is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CuCLARK. If not, see <http://www.gnu.org/licenses/>.


COMPATIBILITY TO CLARK
-----
CuCLARK uses many of the same command parameters as CLARK, but removes some of them and 
introduces new ones (see MANUAL & OPTIONS for details). CuCLARK's only mode of operation 
is comparable to CLARK's "full" mode. For a low-level description of CLARK's main features 
and options please see "README_CLARK.txt". CuCLARK utilizes adjusted versions of CLARK's 
scripts for metagenomic classification (see CLASSIFICATION OF METAGENOMIC SAMPLES).

In contrast to CLARK, CuCLARK uses canonical k-mers in the database and for the classification. 
Because of this the databases created by CLARK cannot be used for CuCLARK. However, CLARK is 
able to use the databases created by CuCLARK, because it queries for each k-mer and its reverse, 
so the canonical k-mer is found as it is one of the two. To avoid confusion of the databases, 
CuCLARK's target definition script (set_targets.sh) creates a different subfolder for CuCLARK's 
database than CLARK's.

NOTE: It is possible to substitute the hash table source files (hashTable_hh.hh, HashTableStorage_hh.hh) 
of CLARK version 1.1.3 with the ones of CuCLARK. This way CLARK can create a database compatible to 
CuCLARK, even on a computer that does not provide any CUDA device.

The scripts provided by CLARK to evaluate its results (`estimate_abundance.sh`, 
`evaluate_density_confidence.sh`, `evaluate_density_gamma.sh`) can also be used for CuCLARK, 
because it uses the same result format. Please refer to CLARK's README.txt for further 
details on these scripts.


SOFTWARE & SYSTEM REQUIREMENTS
-----
1) C++ COMPILER VERSION
The main requirement is a 64-bit operating system (Linux or Mac), and the GNU GCC to
compile version 4.4 or higher. Multi-threading operations are assured by the openmp
libraries. If these libraries are not installed, CuCLARK will run in single-threaded
task. CuCLARK is expected to successfully run on recent Mac and Linux OS.

2) MEMORY

In our tests CuCLARK needed about 146 GB of RAM to build its database of target-specific 31-mers 
from the bacteria genomes of NCBI/RefSeq. The resulting database occupied about 36 GB of 
storage space on the hard drive. To use the database for classification CuCLARK needed about 40 GB of RAM.
CuCLARK-l can build its database and run the classification on computers with just 4 GB of RAM. 
The database for bacteria was less than 600 MB in size on hard disk.

3) CUDA DEVICES

CuCLARK and CuCLARK-l need at least one CUDA-capable device of compute capability 3.0 or higher 
with at least 1 GB of graphic memory to run. CuCLARK can utilize an arbitrary number of such devices 
to speed up the classification. The database is divided into parts depending on the devices' memory sizes. 
If the database does not fit on the devices at once, CuCLARK performs multiple cycles and swaps the 
database parts in and out. For this reason it is advised to provide as much graphic memory as possible. 
In our test we used 1-4 devices with 6 GB of memory each.
CuCLARK-l can fit its database in one part on devices with 1 GB of graphic memory.

3) CUDA VERSION

CuCLARK was tested with CUDA version 7.5. Although we expect CuCLARK to work with lower CUDA versions, 
we do not guarantee it.


INSTALLATION
-----
Copy the whole "CuCLARK" folder to hard disk and execute the installation script (`./install.sh`).
The installer builds binaries (CuCLARK and CuCLARK-l, in the subfolder "exe").

SCRIPTS
-----
In the main folder, you can also notice that several scripts are available.
Especially:
- `set_targets.sh` and `classify_metagenome.sh`: They allow you to classify your metagenomes
against several database(s) (downloaded from NCBI or available "locally" in your disk).
See section "CLASSIFICATION OF METAGENOMIC SAMPLES" for details.

- `download_data.sh`, `download_taxondata.sh` and `make_metadata.sh` are called by `set_targets.sh` to download a specific database and taxonomy tree data from NCBI, and to associate the genomes of the database with the corresponding taxons, respectively. Although it is possible to use these scripts on their own, we recommend to simply use `set_targets.sh` to carry out all necessery steps.

- `download_data.sh` downloads bacteria, viruses or human genomes from NCBI like the original CLARK.

- `download_data.sh` can be replaced with `download_data_newest.sh` or `download_data_release.sh`
to download the newest NCBI RefSeq genomes or the genomes of the latest NCBI RefSeq release. These scripts allow to download any database included in RefSeq like archaea, bacteria, fungi, etc..

- `clean.sh`: This script will delete permanently all data related (generated and 
downloaded) of the database directory defined in set_targets.h.

- `resetCustomDB.sh`: It resets the targets definition with sequences (newly 
added/modified) of the customized database. Any call of this script must be 
followed by a run of set_target.sh.

- `updateTaxonomy.sh`: To download the latest taxonomy data (taxonomy id, accession numbers, etc.) from the NCBI website.



Following is a version of CLARK's usage guide adjusted to CuCLARK's needs.

MANUAL & OPTIONS
-----
CuCLARK offers several options to run the classification. 

A typical command line to run CuCLARK (or CuCLARK-l) looks like:
```
$ ./CuCLARK -k <kmerSize> -T <fileTargets> -t <minFreqTarget> -D <directoryDB/> -O <fileObjects> -R <fileResults> -n <numberofthreads> ...
```
Definitions of parameters:

`-k <kmerSize>`,       	 	k-mer length:	integer, >= 2 and <= 32. 
			 	The default value for this parameter is 31, except for CuCLARK-l (it is 27).

`-T <fileTargets>`,    	 	The targets definition is written in fileTargets: filename.
				This is a two-column file (separated by space, comma or tab), such that, for each line:
				column 1: the filename of a reference sequence
				column 2: the target ID (taxon name, or taxonomy ID, ...) of the reference sequence 
				(See example below).
				This parameter is mandatory.

`-t <minFreqTarget>`,    	minimum of k-mer frequency/occurence for the discriminative k-mers:  integer, >=0.
				The default value is 0. For example, for 1 (or, 2), the program will discard any 
				discriminative k-mer that appears only once (or, less than twice).

`-D <directoryDB/>`,   	 	Directory of the database : pathname.
				This parameter is mandatory.

`-O <fileObjects>`,    	 	file containing objects: filename.
				This parameter is mandatory.

`-P <file1> <file2>`,		Paired-end fastq files: filenames.

`-R <fileResults>`,    	 	file to store results:  filename.
				Results are stored in CSV format in the file <fileResults>.csv (the extension 
				".csv" is automatically added to the filename).
				This parameter is mandatory. 
				
`-n <numberofthreads>`,	 	number of threads: integer >= 1.
				The program runs in parallel for the classification if n > 1.

`-b <numberofbatches>`,		number of batches: integer >= 1.
				The input files are divided into batches. The GPUs process one batch at a time.
				It is recommended to choose a high number of batches to reduce the amount of memory required per batch.

`-d <numberofdevices>`,		number of CUDA devices to use: integer >= 1.
				If unspecified, CuCLARK uses all available CUDA devices.

`-g <iterations>`,		"gap" or number of non-overlapping k-mers to pass when creating the
				database (for CuCLARK-l only). The bigger the gap is, the lower the RAM 
				usage is. The default value is 4, but the user can increase this value if 
				needed to reduce the RAM usage (but this will degrade the sensitivity).

`-s <factor>`,			sampling factor value (for CuCLARK-l only). 
				The higher is this value, the higher are the classification speed
				and the precision. However, our experiments show that the sensitivity
				can be degraded.

`--tsk`,               	 	to request a detailed creation of the database 
				(target specific k-mers files). For each target ID, the program
				will create a file containing all target specific k-mers.
				This option may require a high amount of disk space to complete. 
				This option is not available for CuCLARK-l.

`--extended`			to request the extended output of results. This
				includes hit counts for all targets. When selecting this feature,
				the disk spaced by the results file is significantly increased.

When a filename is required, we recommend absolute paths.


CLASSIFICATION OF METAGENOMIC SAMPLES
-----
We provide several scripts to facilitate the classification in the context of 
metagenomics. CuCLARK can preprocess databases of bacteria, viruses or human (downloaded 
from NCBI) or a customized set of genomes.

First, we present here two scripts, `set_targets.sh` and `./classify_metagenome.sh` that 
work together.


1) Setting and classification

1.1) Step I: Setting targets

After installing CuCLARK (`./install.sh`), the user must create a directory to store all 
reference sequences (bacteria, viruses, human and custom). For all our examples below,
we name this directory path in a generic way `<DIR_DB/>` for clarity. This directory
can be anywhere in your disk(s).

Then, the user must indicate what database(s) to consider for the classification among 
bacteria, viruses, human and/or custom.

For example, only bacteria genomes:
`$ ./set_targets.sh <DIR_DB/> bacteria`

To work with bacteria, viruses and human:
`$ ./set_targets.sh <DIR_DB/> bacteria viruses human`

To classify against a custom database:
The user will need to paste its sequences (fasta files with accession numbers in the 
header, i.e., ">accession.number ..." or ">gi|number|ref|accession.number| ...", 
and one fasta file per reference sequence) in the directory "Custom", inside `<DIR_DB/>`. 
To do so, the user must (1) create the directory "Custom" inside  `<DIR_DB/>` (if it
does not exist yet) (2) copy or move sequences of interest in Custom and (3) run:
`$ ./set_targets.sh <DIR_DB/> custom`


In general case (when the user selects bacteria, viruses and/or human), 
if the directory `<DIR_DB/>` is empty, then the script will download all
the selected database(s) and also data of the taxonomy tree, from NCBI.
Once the sequences are found or downloaded in the directory, the script will build 
the targets for a given taxonomy rank. 

The default taxonomy rank is species. To use a different taxonomy rank, for example, 
genus, the command line is (from the example selecting bacteria, viruses and human):

`$ ./set_targets.sh <DIR_DB/> bacteria viruses human --genus`

In the current version of CuCLARK, the user can choose between six ranks (species to phylum):
`--species` (the default value), `--genus`, `--family`, `--order`, `--class` or `--phylum`.

Indeed, the strength of CuCLARK is to be able to classify quickly and accurately 
metagenomic samples by using only one taxonomy rank. So as a general rule when 
classifying metagenomic reads: 
Consider first the genus or species rank, then if a high proportion of reads 
cannot be classified, reset your targets definition at a higher taxonomy rank 
(e.g., family or phylum).

Once set_targets.sh is finished, the user can proceed to the step II. However, 
if the user wants to modify the selected databases and/or taxonomy rank, then he/she will 
need to run again set_targets.sh with updated parameters before proceeding to step II.

1.2) Step II: Running the classification

The script to run the classification of a metagenomic sample against the database(s) 
previously set in step I is `classify_metagenome.sh`.

For your convenience, this script runs the executable CuCLARK (or CuCLARK-l) and allows you 
to pass few parameters.

For example, say objects to be classified are reads in "sample.fa" (e.g., located in the 
current directory), and results to be stored in "result.csv". A basic command line is:

`$ ./classify_metagenome.sh -O ./sample.fa -R ./result`

As explained in the section "MANUAL & OPTIONS", thanks to identifiers `-O` and `-R`, 
the script will pass the objects file "sample.fa" and results will be stored
in "./result.csv". Objects are classified against the targets and the taxonomy rank
defined by the last execution of `./set_targets.sh`.

IMPORTANT NOTES:
- The targets definition is automatically passed to CLARK in step II. The file has been
 computed by `set_targets.sh`.

- The script `set_targets.sh` assumes that each reference file from bacteria, viruses or custom
database contains an accession number (in the RefSeq records format: 
i.e., ">accession.number ..." or ">gi|number|ref|accession.number| ..." ). 
If a GI number is missing in a file, then this file will not be used for the classification. 

- `set_targets.sh` also maps the accession number found in each reference sequence to its
taxonomy ID based on the latest NCBI taxonomy data. If a mapping cannot be made for a given sequence, 
then it will NOT be counted and excluded from the targets definition.
The total number of excluded files is prompted in the standard output, and all files that have
been excluded are reported in the file "files_excluded.txt" (located in the the specified
database directory (i.e., "./DBD/").
If some files are excluded, then it will probably mean that they have been removed 
for curations for example (visit the RefSeq FAQ webpage).

- You can update your local taxonomy database thanks to the script `updateTaxonomy.sh`
You can use this script before running `set_targets.sh`.

- If the user wants to work with a different customized database (for example, by removing
or adding more sequences of interest in the Custom folder) then the targets definition
must be reset. We made it simple with the script `resetCustomDB.sh`: 
After the sequences in the Custom folder have been updated, just run:
`$ ./resetCustomDB.sh`
Then, run `set_target.sh` with the desired settings.

- The database files (*.ky, *.lb and *.sz) will be created inside some subdirectory of the 
specified database directory in step I (i.e., "./DBD/") by `classify_metagenome.sh`.

- The default values (the k-mer length, the number of threads, etc.) are used 
if not specified by the user, just like indicated in the previous section.

- The script `classify_metagenome.sh` still allows you to pass customized parameters and 
options, similarly to the previous section. `classify_metagenome.sh` follows options defined
in "MANUAL & OPTIONS"(see below some examples). So you can change the k-mer length,
the number of parallel threads, etc.

We present below some examples of customized classification using classify_metagenome.sh.

To use 20-mers (instead of 31-mers):
`$ ./classify_metagenome.sh -O ./sample.fa -R ./result -k 20`

To request 2 threads and 16 batches:
`$ ./classify_metagenome.sh -O ./sample.fa -R ./result -n 2 -b 16`

To limit the number of CUDA devices to use to 1:
`$ ./classify_metagenome.sh -O ./sample.fa -R ./result -d 1`

To request execution with gzipped objects file, and using 8 threads:
`$ ./classify_metagenome.sh -O ./sample.fa.gz -R ./result -n 8 --gzipped`

Another example, for classifying paired-end reads (./sample1.fastq and ./sample2.fastq):
`$ ./classify_metagenome.sh -P ./sample1.fastq ./sample2.fastq -R ./result`


Note:
This script can run CuCLARK-l instead of CuCLARK, for workstations with limited RAM. 
Then, the user can indicate it with the option  `--light`. For example:

`$ ./classify_metagenome.sh -P ./sample1.fastq ./sample2.fastq -R ./result --light`

Note:
Typing only `./classify_metagenome.sh` in the terminal will prompt the help describing 
options and parameters.


RESULTS FORMAT
-----
Results are in CSV format.

The default results format is the following for each line:
`<Object_ID>,<Length of object>,<Gamma>,<first assignment>,<hit count of first>,<second assignment>,<hit count of second>,<confidence score>` where :
* the "Object_ID"        is the tag name indicated in the header (after ">" or "@") for each object.
* Length of object       is the number of bases (A,C,G,T,U or N) the object has.
* Gamma                  is the ratio between the total number of hit found in the object 
			 (against all targets) and the number of k-mers in the object.
* first assignment       is the target ID of the target that obtained the highest hit count
                         (ties are broken arbitrarily).
* hit count of first     is the number of hit count for the first assignment (h1).
* second assignment      is the target ID of the target that obtained the second highest hit 
                         count (ties are broken arbitrarily).
* hit count of second    is the number of hit count for the second assignment (h2).
* confidence score       is the ratio : h1/(h1+h2).


With the option "--extended", the results format is the following for each line:
`<Object_ID>,<hit count in target 1>,...,<hit count in target N>,<Length of object>,<Gamma>,<first assignment>,<hit count of first>,<second assignment>,<hit count of second>,<confidence score>` where :
* the "Object_ID"        is the tag name indicated in the header (after ">" or "@") for each 
			 object, and N is the number of targets. 
* hit count in target i  is the number of k-mers specific to target i that are in the object.
* Length of object       is the number of bases (A,C,G,T,U or N) the object has.
* Gamma                  is the ratio between the total number of hit found in the object 
			 (against all targets) and the number of k-mers in the object.
* first assignment       is the target ID of the target that obtained the highest hit count
                         (ties are broken arbitrarily).
* hit count of first     is the number of hit count for the first assignment (h1).
* second assignment      is the target ID of the target that obtained the second highest hit 
                         count (ties are broken arbitrarily).
* hit count of second    is the number of hit count for the second assignment (h2).
* confidence score       is the ratio : h1/(h1+h2).




