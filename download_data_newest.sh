#!/bin/sh

#
#   cuCLARK, CLARK for CUDA-enabled GPUs.
#   Copyright 2016-2017, Robin Kobus <rkobus@students.uni-mainz.de>
#  
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#   download_data_newest.sh: To download newest RefSeq genomes from NCBI site (for Bacteria, 
#		      	Viruses, etc.).

if [ $# -ne 2 ]; then
	echo "Usage: $0 <Directory for the sequences> <Database: bacteria, viruses, ...> "
	exit
fi

if [ "$2" = "bacteria" ]; then
if [ ! -s $1/.bacteria ]; then
	rm -Rf $1/Bacteria $1/.bacteria.*
	mkdir -m 775 $1/Bacteria
	cd $1/Bacteria/
	echo "Downloading now Bacteria genomes:"
	wget ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
	if [ -s "assembly_summary.txt" ]; then
		awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
		awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
		wget -nc -i ftpfilepaths

		echo "Downloading done. Uncompressing files... "
		gunzip *.gz

		rm -f ftpdirpaths
		rm -f ftpfilepaths
	else
		echo "Error: Couldn't find assembly_summary text file!"
		exit
	fi

	find `pwd` -name '*.fna' > ../.bacteria
	cd ..
	if  [ ! -s .bacteria ]; then
		echo "Error: Failed to download bacteria sequences. "
		exit
	fi
	echo "Bacteria sequences downloaded!"
else
	echo "Bacteria sequences already in $1."
fi
exit

fi

if [ "$2" = "viruses" ]; then
if [ ! -s $1/.viruses ]; then
	rm -Rf $1/Viruses  $1/.viruses.*
	mkdir -m 775 $1/Viruses
	cd $1/Viruses/
	echo "Downloading now Viruses genomes:"
	wget ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt
	if [ -s "assembly_summary.txt" ]; then
		awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
		awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
		wget -nc -i ftpfilepaths

		echo "Downloading done. Uncompressing files... "
		gunzip *.gz

		rm -f ftpdirpaths
		rm -f ftpfilepaths
	else
		echo "Error: Couldn't find assembly_summary text file!"
		exit
	fi

	find `pwd` -name '*.fna'  > ../.viruses
	cd ..
	if  [ ! -s .viruses ]; then
		echo "Error: Failed to download viruses sequences. "
		exit
	fi
	echo "Viruses sequences downloaded!"
else
        echo "Viruses sequences already in $1."
fi
exit
fi

if [ ! -s $1/.$2 ]; then
	rm -Rf $1/$2  $1/.$2.*
	mkdir -m 775 $1/$2
	cd $1/$2/
	echo "Downloading now '$2' genomes:"
	wget ftp://ftp.ncbi.nih.gov/genomes/refseq/$2/assembly_summary.txt
	if [ -s "assembly_summary.txt" ]; then
		awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
		awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
		wget -nc -i ftpfilepaths

		echo "Downloading done. Uncompressing files... "
		gunzip *.gz

		rm -f ftpdirpaths
		rm -f ftpfilepaths
	else
		echo "Error: Couldn't find assembly_summary text file! Are you sure '$2' database exists in RefSeq?"
		exit
	fi

	find `pwd` -name '*.fna'  > ../.$2
	cd ..
	if  [ ! -s .$2 ]; then
		echo "Error: Failed to download '$2' sequences. "
		exit
	fi
	echo "'$2' sequences downloaded!"
else
        echo "'$2' sequences already in $1."
fi
exit
fi

