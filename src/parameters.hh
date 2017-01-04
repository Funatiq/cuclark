/*
   CuCLARK, CLARK for CUDA-enabled GPUs.
   Copyright 2016-2017, Robin Kobus <rkobus@students.uni-mainz.de>
   
   based on CLARK version 1.1.3, CLAssifier based on Reduced K-mers.
   Copyright 2013-2016, Rachid Ounit <rouni001@cs.ucr.edu>


   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @author: Robin Kobus, masters student at Institute of Computer Science, JGU Mainz
 * @project: CuCLARK, Metagenomic Classification with CUDA-enabled GPUs
 * 
 * Changes:
 * Added additional defines.
 */


#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#define VERSION "1.1"

#define SB              4       
#define LHTSIZE 	57777779
#define HTSIZE  	1610612741
#define NBN		1
#define SFACTORMAX 	30

//// new defines
#define MAXHITS		15			// maximum number of targets per object, maximum hits in normal test 10M was 10
#define RESERVED	400000000	// reserverd GPU memory for a batch
#define OBJECTNAMEMAX	40		// maximum length for object names

#define DBPARTSPERDEVICE 3
////

typedef uint64_t      T64;
typedef uint32_t      T32;
typedef uint16_t      T16;

#endif
