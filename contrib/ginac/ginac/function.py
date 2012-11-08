#!/usr/bin/env python
# encoding: utf-8

maxargs = 14
methods = "eval evalf conjugate real_part imag_part derivative power series print".split()

import sys, os, optparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
import yaptu

def seq(template, N, start = 1, separator = ', '):
	return separator.join([ template % { 'n' : n } for n in range(start, N + start) ])

def main():
	parser = optparse.OptionParser()
	parser.add_option("-o", dest="outfile")
	options, args = parser.parse_args()
	cop = yaptu.copier(globals())
	if not options.outfile is None:
		cop.ouf = open(options.outfile, 'wt')
	inp = sys.stdin
	if len(args) >= 1:
		inp = open(args[0])
	cop.copy(inp)
	if inp != sys.stdin:
		inp.close()
	if cop.ouf != sys.stdout:
		cop.ouf.flush()
		cop.ouf.close()
	return 0

if __name__== '__main__':
	ret = main()
	sys.exit(ret)

