#! /bin/sh
#
#   SUMMARY: 
#     USAGE:
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@epfl.ch
#
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

the_test=
the_srcdir=`cat Makefile | grep "^srcdir =" | awk 'BEGIN{FS="="}{ print $2} '`


compile()
{
cd $1 && make -s $1 > /dev/null 2>&1
if test "$?" != "0"; then
 exit 1;
fi
}

run()
{
cd $1 && ./$1 > $1.log 2>&1;
if test "$?" != "0"; then
 exit 1;
fi
}

for ac_option
do
  # If the previous option needs an argument, assign it.
  if test -n "$ac_prev"; then
    eval "$ac_prev=\$ac_option"
    ac_prev=
    continue
  fi

  ac_optarg=`expr "x$ac_option" : 'x[^=]*=\(.*\)'`
  
  # Accept the important Cygnus configure options, so we can diagnose typos.
  
  case $ac_option in
      
      -c)
	  do_compile=1;;
      -e)
	  do_run=1;;
      -d)
	  do_diff=1;;
      -t=*)
	  the_test=$ac_optarg;;
  esac
done

if test "$do_compile" = "1"; then
    #echo "test compilation of $the_test..."
    (compile $the_test)
fi
if test "$do_run" = "1"; then
    #echo "test execution of $the_test..."
    (run $the_test)
fi

if test  "$do_diff" = "1"; then
    
    #echo "test results of $the_test..."
    cmd=`diff $the_test/$the_test.log $the_srcdir/$the_test/main_sample.out`;
    if ! test -z "$cmd"; then
	exit 1;
    fi
fi

exit 0;
