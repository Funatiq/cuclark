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
#   resetDB.sh: To reset the database files created by CLARK 
#

if [ "$1" = "--help" ]; then
	echo "This script erases all database files created with old Custom sequences."
	echo "Please use this script after having updated the Custom folder."
	exit
fi

echo "Are you sure you have updated the Custom directory ? (yes/no)"
read decision

if [ $decision = "yes" ] || [ $decision = "y" ] || [ $decision = "Y" ] || [ $decision = "Yes" ] || [ $decision = "YES" ]; then
echo -n "The program will clean all database files created with the previous data in the Custom directory..."
for DIR in `cat ./.DBDirectory`
do
rm -f $DIR/targets.txt
rm -Rf $DIR/custom*
rm -Rf $DIR/*_custom*
rm -f $DIR/.custom*

echo "done"
echo -n "Resetting the list of custom sequences..."
find $DIR/Custom/ -name '*.f*' > $DIR/.custom
echo "done"
done

else

exit
fi

