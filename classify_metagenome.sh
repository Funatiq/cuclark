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

#
#   classify_metagenome.sh: To classifiy a metagenomic sample against one
#			   or several databases from the list:
#			   bacteria, viruses, human and custom.
#
#   Differences to CLARK:
#   Deprecated some parameters and added new parameters.
#

if [ $# -lt 2 ]; then
echo "Usage: $0 -k <kmerSize> -O <fileObjects> -R <fileResults> -n <numberofthreads> ... "
echo ""
echo "Definitions of parameters (see README for more details):"
echo ""
echo "-k <kmerSize>,       \t k-mer length:\tinteger, >= 2 and <= 32. The default value is 31.\n"
echo "-t <minFreqTarget>,  \t minimum of k-mer frequency in targets:\tinteger, >=0. The default value is 0.\n"
#echo "-o <minFreqtObject>, \t minimum of k-mer frequency in objects:\tinteger, >=0. The default value is 0.\n"
echo "-O <fileObjects>,    \t file containing objects to classify:\t filename.\n"
echo "-P <file1> <file2>,  \t paired-end reads:\t filenames.\n"
echo "-R <fileResults>,    \t file to store results (or corresponding list of results file):\t filename.\n"
#echo "-m <mode>,           \t mode of execution: 0 (full), 1 (default), 2 (express) or 3 (spectrum).\n"
echo "-n <numberofthreads>,\t number of threads:\tinteger >= 1.\n"
echo "-b <numberofbatches>,\t number of batches:\tinteger >= 1.\n";
echo "-d <numberofdevices>,\t number of CUDA devices to use:\tinteger >= 1.\n";
echo "-g <iteration>,      \t gap or number of non-overlapping k-mers to pass for the database creation (for cuCLARK-l only). The default value is 4.\n"
echo "-s <factor>,         \t sampling factor value (for cuCLARK only). \n"
echo "--tsk,               \t to request a detailed creation of the database (target specific k-mers files).\n"
#echo "--ldm,               \t to request the loading of the database by memory mapped-file (in multithreaded mode, multiple parallel threads are requested).\n"
#echo "--kso,               \t to request a preliminary k-spectrum analysis of each object (for mode 3 only).\n"
echo "--extended,	   \t to request an extended output for the results file.\n"
echo "--light		   \t to run the RAM-light variant of cuCLARK, namely cuCLARK-l.\n"
echo "--gzipped            \t to indicate that objects are gzipped.\n"
exit
fi

if [ ! -s .settings ]; then
echo "Please set the targets (cf. the script 'set_targets.sh') before running the classification."
exit
fi


PARAMS=""
UNRQTD=0
FILE1="n"
FILE2="n"
FILE1U=""
FILE2U=""
VARIANT="DEFAULT"
for param in $@
do
	if [ "$param" = "--gzipped" ]; then
                UNRQTD="--gzipped"
                continue
        fi
	if [ "$param" = "--light" ]; then
        	VARIANT="LIGHT"
	fi
done

for line in `cat .settings`
do
	PARAMS="$PARAMS $line"
done

for var in $@
do
	LV=$UNRQTD
	if [ "$var" = "-T" ]; then
		echo "Attempt to overwrite the targets definition: '-T <...>'."
		echo "The targets definition cannot be reset/changed by classify_metagenome.sh." 
		echo "It was previously set by set_targets.sh.\nPlease rerun classify_metagenome.sh without setting a target definition."
		exit
	fi
	if [ "$var" = "-D" ]; then
                echo "Attempt to overwrite the directory of the database: '-D <...>'."
		echo "The database directory cannot be reset/changed in classify_metagenome.sh." 
                echo "It was previously set by set_targets.sh. \nPlease rerun classify_metagenome.sh without setting a database directory."
                exit
        fi
	if [ "$var" = "-O" ]; then
		FILE1="s"
		FILE2="n"
		PARAMS="$PARAMS $var"
		continue
	fi
	if [ "$var" = "-P" ]; then
                FILE1="s"
                FILE2="s"
		PARAMS="$PARAMS $var"
                continue
        fi
	if [ "$FILE1" = "s" ]; then
                FILE1=$var
                if [ "$LV" = "--gzipped" ]; then
			PRE="."
		        FILE1=$PRE$var
			cp $var $FILE1
			FILE1U=`echo $FILE1 | rev | cut -c 4- | rev`
			gunzip -f $FILE1
			PARAMS="$PARAMS $FILE1U"
                else
                        PARAMS="$PARAMS $var"
                fi
                continue
        fi	
	if [ "$FILE2" = "s" ]; then
		if [ "$LV" = "--gzipped" ]; then
                	PRE="."
                        FILE2=$PRE$var
                        cp $var $FILE2
			FILE2U=`echo $FILE2 | rev | cut -c 4- | rev`
			gunzip -f $FILE2
			PARAMS="$PARAMS $FILE2U"
		else
			PARAMS="$PARAMS $var"
		fi
		continue
        fi
	if [ "$var" = "--gzipped" ] || [ "$var" = "--light" ]; then
        	continue
	fi
	PARAMS="$PARAMS $var"
done

if [ ! -f $FILE1U ] || [ ! -f $FILE2U ]; then
echo "Failed to uncompress input objects."
echo "The program must abort."
exit
fi

if [ "$VARIANT" = "LIGHT" ]; then
./exe/cuCLARK-l $PARAMS
else
./exe/cuCLARK $PARAMS
fi

if [ "$UNRQTD" = "--gzipped" ]; then
	if [ "$FILE1U" != "" ]; then
		rm -f $FILE1U
	fi
	if [ "$FILE2U" != "" ]; then
		rm -f $FILE2U
	fi
fi

