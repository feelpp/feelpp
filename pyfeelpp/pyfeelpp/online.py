import pyfeelpp
import pyfeelpp.core
import pyfeelpp.crb
import pyfeelpp.online
import numpy
import pymongo
import sys,time

class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name)
        print('Elapsed: %s' % (time.time() - self.tstart))

def online( plugindir, pluginname, dbid, N=1, load=pyfeelpp.crb.CRBLoad.rb ):
     e=pyfeelpp.core.Environment(sys.argv,pyfeelpp.crb.makeCRBOptions())
     plugin=pyfeelpp.crb.factoryCRBPlugin(pluginname)
     print("dbid:",dbid," name:",plugin.name());
     plugin.loadDBFromId(id=dbid,load=load,root="/home/prudhomm/feel/");
     Dmu=plugin.parameterSpace()
     s=Dmu.sampling()
     s.sampling(N,"random")
     if ( load == pyfeelpp.crb.CRBLoad.all ):
          plugin.initExporter()
     with Timer('plugin.run'):
          r=plugin.run(s.getVector())
     for x in r:
          print("mu=",x.parameter(),"s=",x.output(), " error=",x.errorBound())
          if ( load == pyfeelpp.crb.CRBLoad.all ):
               plugin.exportField("sol-"+str(x.parameter()),x)
     if ( load == pyfeelpp.crb.CRBLoad.all ):
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
     print("members:",pyfeelpp.crb.CRBLoad.__members__)
     print("rb members:",pyfeelpp.crb.CRBLoad.__members__["rb"])

     feelppdb = pymongo.MongoClient()
     crbdb = feelppdb['feelpp']['crbdb']
     
     model=crbdb.find_one({'crbmodel.name':options.name})
     if model == None:
         print("Invalid model name ",options.name)
         exit
     import pprint
     pprint.pprint(model)
     online(pluginname=options.name,plugindir=options.dir, dbid=model['uuid'],N=options.N,load=pyfeelpp.crb.CRBLoad.__members__[options.load])
     
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
    
