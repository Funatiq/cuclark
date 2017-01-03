#!/bin/sh

#
#   cuCLARK, CLARK for CUDA-enabled GPUs.
#   Copyright 2016-2017, Robin Kobus <rkobus@students.uni-mainz.de>
#   
#   based on CLARK version 1.1.3, CLAssifier based on Reduced K-mers.
#   Copyright 2013-2016, Rachid Ounit <rouni001@cs.ucr.edu>
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
#   download_data.sh: To download genomes from NCBI site (for Bacteria, 
#		      Viruses, and Human).

if [ $# -ne 2 ]; then
	echo "Usage: $0 <Directory for the sequences> <Database: bacteria, viruses or human> "
	exit
fi

if [ "$2" = "bacteria" ]; then

if [ ! -s $1/.bacteria ]; then
	rm -Rf $1/Bacteria $1/.bacteria.*
	mkdir -m 775 $1/Bacteria
	cd $1/Bacteria/
	echo "Downloading now Bacteria genomes:"
       	wget ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz
	echo "Downloading done. Uncompressing files... "
	tar -zxf all.fna.tar.gz
	rm all.fna.tar.gz
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
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.ffn.tar.gz
	echo "Downloading done. Uncompressing files... "
	tar -zxf ./all.fna.tar.gz
	tar -zxf ./all.ffn.tar.gz
	rm all.f??.tar.gz
	find `pwd` -name '*.f??'  > ../.viruses
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

if [ "$2" = "human" ]; then
if [ ! -s $1/.human ]; then
	rm -Rf $1/Human  $1/.human.*
	mkdir -m 775 $1/Human
	cd $1/Human/
	echo "Downloading now latest Human genome:"
	n=1

	while [ $n -le 9 ]; do
		wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_0$n/hs_ref_GRC*chr$n.fa.gz
		n=$((n+1))
	done
	n=10
	while [ $n -le 22 ]; do
		wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_$n/hs_ref_GRC*chr$n.fa.gz
		n=$((n+1))
	done

	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_X/hs_ref_GRC*chrX.fa.gz
	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_Y/hs_ref_GRC*chrY.fa.gz
	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_MT/hs_ref_GRC*chrMT.fa.gz
	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_Un/hs_ref_GRC*chrUn.fa.gz
	echo "Downloading done. Uncompressing files... "
	gunzip ./*fa.gz

	find `pwd` -name '*.fa' > ../.human
	cd ../
	if  [ ! -s .human ]; then
                echo "Error: Failed to download human sequences. "
                exit
        fi
	echo "Human sequences downloaded!"
else
        echo "Human sequences already in $1."
fi
exit
fi

echo "Failed to recognize parameter: $2. Please choose between: bacteria, viruses, human."

