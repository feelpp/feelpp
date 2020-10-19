import feelpp
import sys
e=feelpp.Environment(sys.argv)
rd = feelpp.RemoteData("github:{repo:feelpp,path:README.adoc}", worldComm=e.worldCommPtr())

if rd.canDownload():
    d=e.downloadsRepository()
    print("download data in ", d);
    data=rd.download( d )
    print("downloaded data:",data)
