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
#   cleanDB.sh: To erase all databases directories/files generated/downloaded. 
#


if [ ! -s "./.DBDirectory" ]; then
echo "There is no database directory to clean"
exit
fi
DIR=`cat ./.DBDirectory`
echo "Are you sure you want to delete all data in the database directoy: $DIR? (yes/no)"
read decision

if [ $decision = "yes" ] || [ $decision = "y" ] || [ $decision = "Y" ] || [ $decision = "Yes" ] || [ $decision = "YES" ]; then
echo "Cleaning: on-going..."
rm -Rf $DIR
rm -f .dbAddress
rm -f .DBDirectory 
rm -f .settings
echo "Cleaning: done."
elif [ $decision = "no" ] || [ $decision = "n" ] || [ $decision = "N" ] || [ $decision = "No" ] || [ $decision = "NO" ]; then
echo "Cleaning: canceled"
fi
