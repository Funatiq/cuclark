#!/bin/sh

#
#   cuCLARK, CLARK for CUDA-enabled GPUs.
#   Copyright 2016, Robin Kobus <rkobus@students.uni-mainz.de>
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
#   set_targets.sh: To create targets definition of selected databases
#                   (Bacteria, Viruses, Human and Custom).

#   Differences to CLARK:
#   Changed subDB folder name to avoid confusion of CLARK's and cuCLARK's databases.
#

if [ $# -lt 2 ]; then

echo "Usage: $0 <Directory path> <Databases: bacteria, viruses, human or custom>+ <taxonomy rank: --phylum, --class, --order, --family, --genus or --species (default)>"

exit
fi

DBDR=$1
RANK=0

if [ ! -d $DBDR ]; then
	echo "Selected directory not found. The program will create it."
	mkdir -m 775 $DBDR
	if [ ! -d $DBDR ]; then
	echo "Failed to create the directory (please check the name of directory $DBDR and whether it exists). The program will abort."
	exit
	fi
fi
echo $DBDR > .DBDirectory
for var in $@
do
        if [ "$var" = "--species" ]; then
        RANK=0
        break
        fi
        if [ "$var" = "--genus" ]; then
        RANK=1
        break
        fi
        if [ "$var" = "--family" ]; then
        RANK=2
        break
        fi
        if [ "$var" = "--order" ]; then
        RANK=3
        break
        fi
        if [ "$var" = "--class" ]; then
        RANK=4
        break
        fi
        if [ "$var" = "--phylum" ]; then
        RANK=5
        break
        fi
	PREF=`echo $var | cut -c1-2`
        if [ "$PREF" = "--" ]; then
		echo "Failed to recognize this parameter: $var"
		exit 
	fi
done

if [ -f $DBDR/targets.txt ]; then
	rm -f $DBDR/targets.txt
fi

touch $DBDR/targets.txt
rm -f $DBDR/.tmp .settings files_excluded.txt $DBDR/files_excluded.txt
subDB=""
us="_"
for db in $@
do
	if [ "$db" != "$DBDR" ]; then
		PRE=`echo $db | cut -c1-2`
		if [ "$PRE" != "--" ]; then
			echo -n "Collecting metadata of $db... "
			./make_metadata.sh $db $DBDR
			if [ ! -s $DBDR/.$db ]; then
				exit
			fi
			if [ ! -f $DBDR/.taxondata ]; then
				exit
			fi
			echo "done."
			if [ -s $DBDR/.$db.fileToTaxIDs ]; then 
				./exe/getTargetsDef $DBDR/.$db.fileToTaxIDs $RANK >> $DBDR/targets.txt 
				subDB="$subDB$db$us"
				cat files_excluded.txt >> $DBDR/.tmp
				rm files_excluded.txt
			fi
		fi
	fi
done

subDB="$subDB$RANK"
subDB="${subDB}_canonical"
echo "-T $DBDR/targets.txt" > .settings
if [ ! -d $DBDR/$subDB ]; then
	echo "Creating directory to store discriminative k-mers: $DBDR/$subDB"
	mkdir -m 775 $DBDR/$subDB
fi
echo "-D $DBDR/$subDB/" >> .settings
if [ -s $DBDR/.tmp ]; then
	mv $DBDR/.tmp $DBDR/files_excluded.txt
fi
