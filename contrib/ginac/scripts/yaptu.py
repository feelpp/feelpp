#!/usr/bin/env python
"""
Yet Another Python Templating Utility, Version 1.2, by Alex Martelli.
Distributed under PSF license (http://docs.python.org/license.html).
"""

import re
import sys

# utility stuff to avoid tests in the mainline code
class _nevermatch:
    "Polymorphic with a regex that never matches"
    def match(self, line):
        return None
_never = _nevermatch()     # one reusable instance of it suffices
def identity(string, why):
    "A do-nothing-special-to-the-input, just-return-it function"
    return string
def nohandle(string):
    "A do-nothing handler that just re-raises the exception"
    raise

_default_rex = re.compile('@([^@]+)@')
_default_rbe = re.compile('\+\+\+')
_default_ren = re.compile('---')
_default_rco = re.compile('===')

# and now the real thing
class copier:
    "Smart-copier (YAPTU) class"
    def copyblock(self, i=0, last=None):
        "Main copy method: process lines [i,last) of block"
        def repl(match, self=self):
            "return the eval of a found expression, for replacement"
            # uncomment for debug:
	    # print '!!! replacing',match.group(1)
            expr = self.preproc(match.group(1), 'eval')
            try: return str(eval(expr, self.globals, self.locals))
            except: return str(self.handle(expr))
        block = self.locals['_bl']
        if last is None: last = len(block)
        while i<last:
            line = block[i]
            match = self.restat.match(line)
            if match:   # a statement starts "here" (at line block[i])
                # i is the last line to _not_ process
                stat = match.string[match.end(0):].strip()
                j=i+1   # look for 'finish' from here onwards
                nest=1  # count nesting levels of statements
                while j<last:
                    line = block[j]
                    # first look for nested statements or 'finish' lines
                    if self.restend.match(line):    # found a statement-end
                        nest = nest - 1     # update (decrease) nesting
                        if nest==0: break   # j is first line to _not_ process
                    elif self.restat.match(line):   # found a nested statement
                        nest = nest + 1     # update (increase) nesting
                    elif nest==1:   # look for continuation only at this nesting
                        match = self.recont.match(line)
                        if match:                   # found a contin.-statement
                            nestat = match.string[match.end(0):].strip()
                            stat = '%s _cb(%s,%s)\n%s' % (stat,i+1,j,nestat)
                            i=j     # again, i is the last line to _not_ process
                    j=j+1
                stat = self.preproc(stat, 'exec')
                stat = '%s _cb(%s,%s)' % (stat,i+1,j)
                # for debugging, uncomment...:
		# print "-> Executing: {"+stat+"}"
                exec stat in self.globals,self.locals
                i=j+1
            else:       # normal line, just copy with substitution
                self.ouf.write(self.regex.sub(repl,line))
                i=i+1
    def __init__(self, dict={}, regex=_default_rex,
            restat=_default_rbe, restend=_default_ren, recont=_default_rco,
            preproc=identity, handle=nohandle, ouf=sys.stdout):
        "Initialize self's attributes"
        self.regex   = regex
        self.globals = dict
        self.locals  = { '_cb':self.copyblock }
        self.restat  = restat
        self.restend = restend
        self.recont  = recont
        self.preproc = preproc
        self.handle  = handle
        self.ouf     = ouf
    def copy(self, inf=sys.stdin, block=None):
        "Entry point: copy-with-processing a file, or a block of lines"
        if block is None: block = inf.readlines()
        self.locals['_bl'] = block
        self.copyblock()


