#!/usr/bin/python

from xmlParser import *
from xml.dom import minidom

################################################################################
#				Main program
################################################################################
params,outputs = parse_xml_response("essai_xml_response2.xml")
for p in params:
	print p.getName()
print "==================="
print outputs

data = parse_xml_result("essai_xml2.xml", ["program","dim","order","beta","nu","h"], ["laplacian","1","1","0.1","0.2","0.5"], "norm_L2")
print data