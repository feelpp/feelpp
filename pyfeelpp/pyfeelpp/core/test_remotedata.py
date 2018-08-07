import pyfeelpp.core as core
import sys
e=core.Environment(sys.argv)
rd = core.RemoteData("github:{repo:feelpp,path:README.adoc}", worldComm=e.worldComm())
if rd.canDownload():
    d=e.downloadsRepository()
    print("download data in ", d);
    rd.download( d )
