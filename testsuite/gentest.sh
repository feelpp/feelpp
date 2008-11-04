#! /bin/sh
#
#   SUMMARY: generate test
#     USAGE: gentest
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@epfl.ch
#
# DESCRIPTION:
# ============
# Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# DESCRIP-END.
#

for i in $1
do
  if test -d $i; then
      cp /dev/null  $i/testsuite.at
      keywords=`echo $i | awk 'BEGIN{FS="_"}{for(i=2;i<=NF;++i) print $i;}'|xargs echo -n`
      echo -n "generating $i/testsuite.at with keywords=$keywords..."
      echo "AT_SETUP([$i])" >> $i/testsuite.at
      echo "AT_KEYWORDS([$keywords])" >> $i/testsuite.at
      if test -f $i/data; then
	  top_builddir=`cat Makefile.in| grep top_builddir | head -1 | sed "s/top_builddir = //g"`
	  
	  meshfile=`cat $i/data | grep mesh_file | head -1 | sed "s/ //g" | awk 'BEGIN{FS="="}{print $2}'|awk 'BEGIN{FS="#"}{print $1}'`
	  meshdir=`cat $i/data | grep mesh_dir | head -1 | sed "s/ //g" |sed "s/\.\.\///g" | awk 'BEGIN{FS="="}{print $2}'|awk 'BEGIN{FS="#"}{print $1}'`
	  datafile=`cat $i/data | sed "s/mesh_dir[ \t]*=.*$/mesh_dir = .\//g"`

	  cat >> $i/testsuite.at <<EOF
AT_DATA([data],[[$datafile
]])
EOF
      fi
      echo "mesh: $meshdir/$meshfile"
      if test -z $meshfile -a -z $meshdir; then

	  echo "AT_CHECK([\$top_builddir/testsuite/$i/$i],[0],[ignore],[ignore])" >> $i/testsuite.at

      else
	  echo "AT_CHECK([ln -sf \$top_builddir/testsuite/$meshdir/$meshfile &&  \$top_builddir/testsuite/$i/$i],[0],[ignore],[ignore])" >> $i/testsuite.at
      fi
      echo "AT_CLEANUP" >> $i/testsuite.at
      echo "done."
  fi
done
