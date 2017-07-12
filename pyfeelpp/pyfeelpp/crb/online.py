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


def main(args):
     from optparse import OptionParser
     parser = OptionParser()
     parser.add_option('-N', '--Num', dest='N', help='number of parameters to evaluate the plugin', type=int, default=1)
     parser.add_option('-d', '--dir', dest='dir', help='directory location of the plugins', default="/usr/local/bin")
     parser.add_option('-n', '--name', dest='name', help='name of the plugin',type="string")
     parser.add_option('-i', '--id', dest='dbid', help='DB id to be used', type="string")
     parser.add_option('-l', '--load', dest='load', help='type of data to be loadedoaded', type="string",default="rb")
     (options, args) = parser.parse_args()
     print crb.CRBLoad.__members__
     print crb.CRBLoad.__members__["rb"]

     online(pluginname=options.name,plugindir=options.dir, dbid=options.dbid,N=options.N,load=crb.CRBLoad.__members__[options.load])
     
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
    
