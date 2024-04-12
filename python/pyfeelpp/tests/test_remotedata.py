import feelpp.core as fppc
import sys
import pytest

#@pytest.mark.mpi
def test_remotedata(init_feelpp):
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/core/test_remotedata")
    rd = fppc.RemoteData("github:{repo:feelpp,path:README.adoc}", worldComm=fppc.Environment.worldCommPtr())

    if rd.canDownload():
        d=fppc.Environment.downloadsRepository()
        print("download data in ", d);
        data=rd.download( d )
        print("downloaded data:",data)
