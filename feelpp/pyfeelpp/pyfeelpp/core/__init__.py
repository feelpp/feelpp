import sys,time
from _core import *

def download(data,worldComm):
    """Download remote data file"""
    rd = RemoteData(data,worldComm)
    if rd.canDownload():
        d=Environment.downloadsRepository()
        return rd.download( d )
    else:
        raise RuntimeError("Remote Data " + data + " cannot be downloaded")
