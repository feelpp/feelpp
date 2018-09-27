#!/usr/bin/env python
# encoding: utf-8
# Convert help for ginsh operators from man page to C source
import sys, re, optparse

rxStart = re.compile('^.*GINSH_OP_HELP_START$')
rxEnd = re.compile('^.*GINSH_OP_HELP_END$')
fcnRx = re.compile('^[.]B\s+')
minusRx = re.compile('\\\\[-]')
codeFmt = 'insert_help("operators", "%s\\t" "%s");\n'

def extractHelp(inp, out):
	sym, synopsis = None, None
	seenStart = False
	for l in inp:
		l = l.strip()
		if not seenStart:
			if rxStart.match(l):
				seenStart = True
			continue
		if rxEnd.match(l):
			if sym is not None:
				out.write(codeFmt % ( sym, synopsis ))
			break
		if fcnRx.match(l):
			l = fcnRx.sub('', l)
			sym = minusRx.sub('-', l)
		elif l.lower() == '.tp':
			out.write(codeFmt % ( sym, synopsis ))
			sym, synopsis = None, None
		else:
			synopsis = l

def main():
	op = optparse.OptionParser()
	op.add_option('-o', dest = 'outfile')
	options, args = op.parse_args()
	outfile = sys.stdout
	infile = sys.stdin
	if not options.outfile is None:
		outfile = open(options.outfile, 'wt')
	if len(args) >= 1:
		infile = open(args[0])
	extractHelp(infile, outfile)
	if infile != sys.stdin:
		infile.close()
	outfile.flush()
	if outfile != sys.stdout:
		outfile.close()

if __name__ == "__main__":
	main()
	sys.exit(0)

