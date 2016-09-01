#!/bin/sh

# 
#   CLARK, CLAssifier based on Reduced K-mers.
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
#   Copyright 2013-2015, Rachid Ounit <rouni001@cs.ucr.edu>
#   make_metadata.sh: To create meta-data for the selected database (Bacteria,
# 		      Viruses, Human or Custom). 
#

if [ $# -lt 2 ]; then

echo "Usage: $0 <Database name: bacteria, viruses or human> <Directory name for the database>"

exit
fi

DB=$1
DBDR=$2
TAXDR="taxonomy"

if [ ! -d $DBDR ]; then
	echo "Selected directory not found. The program will create it."
	mkdir -m 775 $DBDR
	mkdir -m 775 $DBDR/Custom

	if [ ! -d $DBDR ]; then
		echo "Failed to find the directory (please check the name of directory $DBDR and whether it exists). The program will abort."
		exit
	fi
else
	if [ ! -d $DBDR/Custom ]; then
		mkdir -m 775 $DBDR/Custom
	fi
fi

if [ ! -d $DBDR/$TAXDR ]; then
	echo "Taxonomy data missing. The program will download data to $DBDR/$TAXDR."
	mkdir -m 775 $DBDR/$TAXDR
	./download_taxondata.sh $DBDR/$TAXDR
fi
if [ ! -f $DBDR/.taxondata ]; then
        echo "Failed to find taxonomy files. The program will try to download them..."
        ./download_taxondata.sh $DBDR/$TAXDR
        if [ ! -f $DBDR/.taxondata ]; then
                echo "Failed to find taxonomy files."
                echo "The program must abort."
                exit
        fi
fi 

if [ "$DB" != "custom" ] && [ "$DB" != "bacteria" ] &&  [ "$DB" != "viruses" ] && [ "$DB" != "human" ]; then
	echo "Failed to recognize the database: '$DB'. "
	echo "-> Supported databases are 'bacteria', 'viruses' and 'human' from NCBI (Please check any typo)."
	echo "-> Eventually, the program can build a database from your own sequences (*.fna files) in your disk. "
	echo "If you want CLARK to use a customized database then please do the following directions: "
	echo "1) Move your sequences (each fasta file defined with a GI number) to $DBDR/Custom/"
	echo "2) Run again this command with the option 'custom' "
	exit
fi

if [ ! -d $DBDR ]; then
	echo "The directory $DBDR does not exit. The program will create it."
	mkdir -m 775 $DBDR
	if [ ! -d $DBDR ]; then
	echo  "Failed to create the directory $DBDR. The program must abort." 
	fi
fi

if [ ! -s $DBDR/.$DB ]; then
	if [ "$DB" != "custom" ]; then
		echo "Sequences for $DB not found. The program will download them."
		./download_data.sh $DBDR $DB
	else
		ls $DBDR/Custom/ > $DBDR/.$DB
		if [ ! -s $DBDR/.$DB ]; then
			echo "The database directory 'Custom' is empty."
			exit
		fi
		find $DBDR/Custom/ -name '*.f*' > $DBDR/.$DB
	fi
fi

if [ ! -f ./exe/getfilesToTaxNodes ] || [ ! -f ./exe/getGInTaxID ]; then
	echo "Something wrong occured (source code may be missing or unusable). The program must abort"
	exit
fi

if [ ! -s $DBDR/.$DB ]; then
	echo "Failed to find the downloaded $DB sequences."
	echo "The program must abort."
	exit
fi

if [ $DB = "human" ]; then
	
	if [ ! -s $DBDR/.$DB.fileToTaxIDs ]; then
		for file in `cat $DBDR/.$DB`
		do
		echo "$file 9606 9605 9604 9443 40674 7711" >> $DBDR/.$DB.fileToTaxIDs
		done
	fi
	exit
fi

if [ ! -s $DBDR/.$DB.fileToGInTaxID ] ; then
	./exe/getGInTaxID $DBDR/.$DB $DBDR/$TAXDR/gi_taxid_nucl.dmp > $DBDR/.$DB.fileToGInTaxID
fi
if [ ! -s $DBDR/.$DB.fileToTaxIDs ]; then
	echo "$DB: Retrieving taxonomy nodes for each sequence based on taxon ID..."
	./exe/getfilesToTaxNodes $DBDR/$TAXDR/nodes.dmp $DBDR/.$DB.fileToGInTaxID > $DBDR/.$DB.fileToTaxIDs
fi
exit

