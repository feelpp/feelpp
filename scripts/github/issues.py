#!/usr/bin/pyhton

import json

import pycurl
from StringIO import StringIO

buffer = StringIO()
c = pycurl.Curl()
c.setopt(c.URL, 'https://api.github.com/repos/feelpp/feelpp/issues?milestone=14&state=closed')
c.setopt(c.WRITEDATA, buffer)
c.perform()
c.close()

body = buffer.getvalue()
# Body is a string in some encoding.
# In Python 2, we can print it without knowing what the encoding is.
print(body)

#system "curl https://api.github.com/repos/feelpp/feelpp/issues\?milestone\=14\&state\=closed > issues.json"

#with open("issues.json") as of:
data = json.loads(body)

for issue in data:
    t = issue['title']
    n = issue['number']
    url = issue['html_url']
    print "   - %s [Issue <a href=\"%s\">%s</a>]" % (t, url, n)
