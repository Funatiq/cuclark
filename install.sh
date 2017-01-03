#! /bin/sh

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
#   install.sh: Install cuCLARK and cuCLARK-l,
#		also programs for target definition
#		in folder ./exe/

echo "Installing cuCLARK, cuCLARK-l and programs for target definition in folder ./exe/"
echo ""
make
