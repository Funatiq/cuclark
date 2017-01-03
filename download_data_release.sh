#!/bin/sh

#
#   cuCLARK, CLARK for CUDA-enabled GPUs.
#   Copyright 2016-2017, Robin Kobus <rkobus@students.uni-mainz.de>
#
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
#   download_data_release.sh: To download genomes of latest RefSeq release from NCBI site 
#			(for Bacteria, Viruses, etc.).

if [ $# -ne 2 ]; then
	echo "Usage: $0 <Directory for the sequences> <RefSeq Database: bacteria, viral, ...> "
	exit
fi

if [ "$2" = "bacteria" ]; then
if [ ! -s $1/.bacteria ]; then
	rm -Rf $1/Bacteria $1/.bacteria.*
	mkdir -m 775 $1/Bacteria
	cd $1/Bacteria/
	wget ftp://ftp.ncbi.nih.gov/refseq/release/RELEASE_NUMBER
	relnum=$(cat RELEASE_NUMBER)
	echo "RefSeq release $relnum found."
	echo "Downloading now Bacteria genomes:"
	wget ftp://ftp.ncbi.nih.gov/refseq/release/bacteria/bacteria.*.genomic.fna.gz
	echo "Downloading done. Uncompressing files... "
	gunzip *.gz
	
	echo "Creating single file for each genome... "
	sed -i 's/\(gi|[0-9]*|ref|\)\([[:graph:]]*\)|/\2/' bacteria.*.genomic.fna
	awk '/^>/ {close(file); file=sprintf("%s.fna",substr($1,2,length($1)-1)); print > file; next;} { print >> file; }' bacteria.*.genomic.fna
	rm -f bacteria.*.genomic.fna
	
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
	wget ftp://ftp.ncbi.nih.gov/refseq/release/RELEASE_NUMBER
	relnum=$(cat RELEASE_NUMBER)
	echo "RefSeq release $relnum found."
	echo "Downloading now Viruses genomes:"
	wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz
	echo "Downloading done. Uncompressing files... "
	gunzip *.gz

	echo "Creating single file for each genome... "
	sed -i 's/\(gi|[0-9]*|ref|\)\([[:graph:]]*\)|/\2/' viral.*.genomic.fna
	awk '/^>/ {close(file); file=sprintf("%s.fna",substr($1,2,length($1)-1)); print > file; next;} { print >> file; }' viral.*.genomic.fna
	rm -f viral.*.genomic.fna

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
	wget ftp://ftp.ncbi.nih.gov/refseq/release/RELEASE_NUMBER
	relnum=$(cat RELEASE_NUMBER)
	echo "RefSeq release $relnum found."
	echo "Downloading now '$2' genomes:"
	wget ftp://ftp.ncbi.nih.gov/refseq/release/$2/$2.*.genomic.fna.gz

	if [ ! -e $2.*.genomic.fna.gz ]; then
		echo "Error: Failed to download '$2' sequences. Are you sure '$2' database exists in RefSeq?"
		cd ..
		rm -rf $2/
		exit
	fi

	echo "Downloading done. Uncompressing files... "
	gunzip *.gz

	echo "Creating single file for each genome... "
	sed -i 's/\(gi|[0-9]*|ref|\)\([[:graph:]]*\)|/\2/' $2.*.genomic.fna
	awk '/^>/ {close(file); file=sprintf("%s.fna",substr($1,2,length($1)-1)); print > file; next;} { print >> file; }' $2.*.genomic.fna
	rm -f $2.*.genomic.fna

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

