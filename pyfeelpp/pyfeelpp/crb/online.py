from pyfeelpp import *
import numpy
import sys,time

class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart)

def online( plugindir, pluginname, dbid, N=1, load=crb.CRBLoad.rb ):
     e=core.Environment(sys.argv,crb.makeCRBOptions())
     plugin=crb.factoryCRBPlugin(plugindir,pluginname)
     plugin.loadDBFromId(id=dbid,load=load);
     Dmu=plugin.parameterSpace()
     s=Dmu.sampling()
     s.sampling(N,"random")
     if ( load == crb.CRBLoad.all ):
          plugin.initExporter()
     with Timer('plugin.run'):
          r=plugin.run(s.getVector())
     for x in r:
          print "mu=",x.parameter(),"s=",x.output(), " error=",x.errorBound()
          if ( load == crb.CRBLoad.all ):
               plugin.exportField("sol-"+str(x.parameter()),x)
     if ( load == crb.CRBLoad.all ):
          plugin.saveExporter()      
     return r,plugin
