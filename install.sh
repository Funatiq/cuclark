#! /bin/sh

#
#   CuCLARK, CLARK for CUDA-enabled GPUs.
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
#   install.sh: Install CuCLARK and CuCLARK-l,
#		also programs for target definition.
# 

## Checking if OPENMP libraries are installed
echo |cpp -fopenmp -dM |grep -i open > .tmp
NB=`wc -l < .tmp`

if [ ! -d ./exe/ ]; then
	mkdir ./exe/
fi

# Install programs for target definition
g++ -o ./exe/getTargetsDef ./src_original/getTargetsDef.cc ./src_original/file.cc -O3
g++ -o ./exe/getGInTaxID ./src_original/getGInTaxID.cc ./src_original/file.cc -O3
g++ -o ./exe/getfilesToTaxNodes ./src_original/getfilesToTaxNodes.cc ./src_original/file.cc -O3

# Install CuCLARK and CuCLARK-l
rm -Rf ./.dCuCLARK/
mkdir ./.dCuCLARK/
cp ./src_original/*.hh ./.dCuCLARK/
cp ./src_original/analyser.cc  ./src_original/file.cc  ./src_original/kmersConversion.cc ./.dCuCLARK/
cp ./src_modified_and_new/*.hh ./.dCuCLARK/
cp ./src_modified_and_new/main.cc ./.dCuCLARK/
cp ./src_modified_and_new/CuClarkDB.* ./.dCuCLARK/

if [ $NB -eq 1 ]; then
	# OPENMP is likely supported
	# Building CuCLARK
	nvcc -Xcompiler -fopenmp -arch=sm_30 -o ./.dCuCLARK/CuCLARK -O3 ./.dCuCLARK/*.cc ./.dCuCLARK/*.cu

    cp ./src_modified_and_new/parameters_light_hh ./.dCuCLARK/parameters.hh
	# Building CuCLARK-l (light version)       
	nvcc -Xcompiler -fopenmp -arch=sm_30 -o ./.dCuCLARK/CuCLARK-l -O3 ./.dCuCLARK/*.cc ./.dCuCLARK/*.cu
else
	# Building CuCLARK
	nvcc -o ./.dCuCLARK/CuCLARK -O3 ./.dCuCLARK/*.cc ./.dCuCLARK/*.cu

    cp ./src_modified_and_new/parameters_light_hh ./.dCuCLARK/parameters.hh
    # Building CuCLARK-l (light version)
	nvcc -o ./.dCuCLARK/CuCLARK-l -O3 ./.dCuCLARK/*.cc ./.dCuCLARK/*.cu
fi
mv ./.dCuCLARK/CuCLARK ./exe/
mv ./.dCuCLARK/CuCLARK-l ./exe/
rm -Rf ./.dCuCLARK .tmp
