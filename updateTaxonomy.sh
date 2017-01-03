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
#   Copyright 2013-2016, Rachid Ounit <rouni001@cs.ucr.edu>
#   updateTaxonomy.sh: To download latest files of taxonomy tree data from NCBI site. 
#


for DIR in `cat ./.DBDirectory`
do
cd $DIR/taxonomy/

# Download taxonomy tree info (GI <-> TaxID, and TaxID: info)
# Download taxonomy tree info (AccessionID <-> TaxID, and TaxID: nodes, merged, names)

echo "Downloading... "
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extrat downloaded data
if [ -s nucl_gb.accession2taxid.gz ] && [ -s taxdump.tar.gz ] && [ -s nucl_wgs.accession2taxid.gz ] ; then
        echo "Uncompressing files... "
        #gunzip gi_taxid_nucl.dmp.gz
        gunzip nucl_wgs.accession2taxid.gz
        gunzip nucl_gb.accession2taxid.gz
        tar -zxf taxdump.tar.gz
        if [ -s nucl_gb.accession2taxid ] && [ -s nodes.dmp ] && [ -s nucl_wgs.accession2taxid ]; then
                cat nucl_gb.accession2taxid > ./nucl_accss
                cat nucl_wgs.accession2taxid >> ./nucl_accss
                touch ../.taxondata
                exit
        else
                echo "Failed to uncompress taxonomy data."
        fi
else
        echo "Failed to download taxonomy data!"
        exit
fi

done
