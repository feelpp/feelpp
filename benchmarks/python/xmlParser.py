# -*- mode: python -*-
#
#  This file is part of the Feel library
#
#  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
#        Date: 2009-04-07
#
#   Copyright (C) 2009
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Universite Joseph Fourier (Grenoble I)
#
# \file codeCalculParaview.py
# \author Florent Vielfaure <florent.vielfaure@gmail.com>
# \date 2009-04-07
#

## \package xmlParser
#  This package permits to read and write xml files to communicate
#  with the c++ codes

from xml.dom import minidom
from numpy import *
from util import *


##	parse_xml_response
#	Parses an XML response file to extract the c++ code capabilities
#
#	@param filename the name of the XML file
#
#	@return pname the name of the c++ code
#	@return params the list of the parameters
#	@return outputs the list of the outputs
def parse_xml_response(filename):
	print("[xmlParser] : opening : "+filename+"...")
	xmldoc = minidom.parse(filename)
	response = xmldoc.firstChild

	pname=getNodeValue(find_node_by_tagname(response,"program_name")[0])

	params=[]
	outputs=[]
	for node in find_node_by_tagname(response,"param"):
		name,attrNames,attrValues=parseNode(node)
		value = getNodeValue(node)
		params.append(Parameter(attrNames,attrValues,value))

	for node in find_node_by_tagname(response,"output"):
		name,attrNames,attrValues=parseNode(node)
		depnam=[]
		depval=[]
		for childnode in node.getElementsByTagName("depend"):
			depnam.append(childnode.getAttribute("value"))
			depval.append(getNodeValue(childnode))
		outputs.append(Output(attrNames,attrValues,depnam,depval))

	return pname,params,outputs

##	parse_xml_result
#	Parses an XML result file to extract the values and the theoritical error
#
#	@param filename the name of the XML file
#	@param name the name of the c++ code
#	@param params the list of the parameters
#	@param values the values associated to the parameters
#	@param output the selected output
#
#	@return data the value for the parameters values and the selected output
def parse_xml_result(filename, name, params, values, output):
	print("[xmlParser] : opening : "+filename+"...")
	xmldoc = minidom.parse(filename)
	tests = xmldoc.firstChild

	node=find_node_by_name(tests,"program","value",name)
	newNode=node
	for ind in range(0,len(params)):
		newNode = find_node_by_name(node,params[ind],"value",values[ind])
		node=newNode
	out=find_node_by_tagname(node,output)[0]

	return getNodeValue(out)

##	find_node_by_name
#	Finds childNodes in the XML tree having a specific tagname, attrname and
#	attrvalue
#
#	@param parentnode the name of the parent node
#	@param tagname the name of the tags which have to be found
#	@param attrname the name of the attibute tag
#	@param attrvalue the name of the attibute
#
#	@return childnode the node found
def find_node_by_name(parentnode,tagname,attrname,attrvalue):
	for childnode in parentnode.getElementsByTagName(tagname):
		if childnode.getAttribute(attrname) == attrvalue:
			return childnode

##	find_node_by_tagname
#	Finds childNodes in the XML tree having a specific tagname
#
#	@param parentnode the name of the parent node
#	@param tagname the name of the tags which have to be found
#
#	@return childnode the node found
def find_node_by_tagname(parentnode,tagname):
	result=[]
	for childnode in parentnode.getElementsByTagName(tagname):
		result.append(childnode)
	return result

##	parseNode
#	return the attributes (names and values) of a node and his name
#
#	@param node the name of the node
#
#	@return the node
#   @return attrNames the names of the attributes
#   @return attrVals the values of the attributes
def parseNode(node):
	name=node.tagName
	attrNames=[]
	attrVals=[]
	for ind in range(0,node.attributes.length):
		attrNames.append(node.attributes.item(ind).nodeName)
		attrVals.append(node.getAttribute(attrNames[ind]))
	return name,attrNames,attrVals

##	find_node_by_tagname
#	return data contain in the node
#
#	@param node the node
#
#	@return node.firstChild.data
def getNodeValue(node):
	return node.firstChild.data

# Local Variables:
# indent-tabs-mode: t
# End:
