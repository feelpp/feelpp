#!/usr/bin/env python
# encoding: utf-8
# TeX Live 2012 seems to dislike files produces by doxygen (1.8.x.y)
# In particular, makeindex(1) program creates invalid index entries like
# \hyperpage{NNN_}
# (note the trailing underscore in the page number). This breaks automatic
# builds and is very annoying. Hence this script. It replaces (broken)
# \hyperpage{NNN_} with \hyperpage{NNN}.
# Note: this is an ugly work around, a proper fix is welcome.
import sys, os, re

def fixupind(fname):
	""" Fix \\hyperpage{NNN_} entries in the ind file @var{fname} """
	tmpout = fname + '.tmp' 
	inp = open(fname)
	out = open(tmpout, 'wt')
	rx = re.compile('(hyperpage)[{]([0-9]+)[_][}]')
	for line in inp:
		out.write(re.sub(rx, '\\1{\\2}', line))
	out.flush()
	out.close()
	inp.close()
	os.rename(tmpout, fname)

if __name__ == '__main__':
	if len(sys.argv) <= 1:
		sys.exit(1)
	fixupind(sys.argv[1])
	sys.exit(0)

