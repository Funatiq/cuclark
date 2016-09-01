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
#   download_taxondata.sh: To download taxonomy tree data from NCBI site. 
#

if [ $# -ne 1 ]; then

echo "Usage: $0 <Directory: directory to store taxonomy data> "
echo "Note: if the chosen directory is not empty, then its content will be erased."
exit

fi

rm -Rf $1
mkdir -m 775 $1

cd $1

# Download taxonomy tree info (GI <-> TaxID, and TaxID: info)
echo "Downloading... "
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extrat downloaded data
if [ -s gi_taxid_nucl.dmp.gz ] && [ -s taxdump.tar.gz ]; then
	echo "Uncompressing files... "
	gunzip gi_taxid_nucl.dmp.gz
	tar -zxf taxdump.tar.gz
	if [ -s gi_taxid_nucl.dmp ] && [ -s nodes.dmp ]; then
		touch ../.taxondata
		exit
	fi
else
	echo "Failed to download taxonomy data!"
	exit
fi

