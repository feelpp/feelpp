#!/usr/bin/env python
# encoding: utf-8
# Convert help for ginsh functions from man page to C source
import sys, re, optparse

rxStart = re.compile('^.*GINSH_FCN_HELP_START$')
rxEnd = re.compile('^.*GINSH_FCN_HELP_END$')
fcnRx = re.compile('^[.]BI\s+')
hlpRx = re.compile('\\\\[-]')
codeFmt = 'insert_help("%s",\n"%s"\n);\n'

def parseProto(pStr):
	pStr = fcnRx.sub('', pStr)
	pStr = re.sub('\s', '', pStr)
	pStr = re.sub('"', '', pStr)
	pStr = re.sub(',', ', ', pStr)
	pStr = re.sub('\\[', ' [', pStr)
	name = pStr.split('(')[0]
	return name, pStr

def extractHelp(inp, out):
	name, proto, synopsis = None, None, None
	seenStart = False
	for l in inp:
		l = l.strip()
		if not seenStart:
			if rxStart.match(l):
				seenStart = True
			continue
		if rxEnd.match(l):
			break
		if fcnRx.match(l):
			name, proto = parseProto(l)
		elif hlpRx.match(l):
			l = hlpRx.sub('', l).strip()
			l = re.sub('"', "'", l)
			synopsis = '%s"\n" - %s' % ( proto, l )
		elif l.lower() == '.br':
			synopsis = synopsis or proto
			out.write(codeFmt % ( name, synopsis ))
			name, proto, synopsis = None, None, None

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

if __name__ == '__main__':
	main()
	sys.exit(0)

