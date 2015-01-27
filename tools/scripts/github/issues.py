#!/usr/bin/python

# change the milestone/release number to retrieve issues
# e.g. 14 -> v0.99.0 (next milestone will be 15)
# milestone are in inscreasing order
import json

#system "curl https://api.github.com/repos/feelpp/feelpp/issues\?milestone\=14\&state\=closed > issues.json"
#use pycurl to retrieve issues
import pycurl
from StringIO import StringIO

buffer = StringIO()
c = pycurl.Curl()
c.setopt(c.URL, 'https://api.github.com/repos/feelpp/feelpp/issues?milestone=15&state=closed')
c.setopt(c.WRITEDATA, buffer)
c.perform()
c.close()

body = buffer.getvalue()
#print(body)


data = json.loads(body)
for issue in data:
    t = issue['title']
    n = issue['number']
    url = issue['html_url']
    print "   - %s [Issue <a href=\"%s\">%s</a>]" % (t, url, n)
