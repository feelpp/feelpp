import feelpp
import sys
import pytest

#@pytest.mark.mpi
def test_remotedata(init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/core/test_remotedata")
    rd = feelpp.RemoteData("github:{repo:feelpp,path:README.adoc}", worldComm=feelpp.Environment.worldCommPtr())

    if rd.canDownload():
        d=feelpp.Environment.downloadsRepository()
        print("download data in ", d);
        data=rd.download( d )
        print("downloaded data:",data)
