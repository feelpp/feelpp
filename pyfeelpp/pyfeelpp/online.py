import pyfeelpp
import pyfeelpp.core
import pyfeelpp.crb
import pyfeelpp.online
import numpy
import pymongo
import pprint
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

pyfeelpp.core.Environment(sys.argv,pyfeelpp.crb.makeCRBOptions())        
def initOnline():
    global feelppenv
    feelppenv=pyfeelpp.core.Environment(sys.argv,pyfeelpp.crb.makeCRBOptions())
    
def listPlugins( plugindir=pyfeelpp.core.Info.libdir() ):
    feelppdb = pymongo.MongoClient()
    crbdb = feelppdb['feelpp']['crbdb']

    names=[]
    for m in crbdb.find():
        names.append(m['crbmodel']['name'])
    return names

def loadAllPlugins( plugindir=pyfeelpp.core.Info.libdir() ):
    feelppdb = pymongo.MongoClient()
    crbdb = feelppdb['feelpp']['crbdb']

    names=[]
    plugins=dict()
    for m in crbdb.find():
        pluginname=m['crbmodel']['name']
        names.append(pluginname)
        plugins[pluginname] = loadPlugin(pluginname)
    return plugins

def loadPlugin( pluginname, plugindir=pyfeelpp.core.Info.libdir(), dbid=None, load=pyfeelpp.crb.CRBLoad.rb ):
    if dbid == None:
        feelppdb = pymongo.MongoClient()
        crbdb = feelppdb['feelpp']['crbdb']
     
        model=crbdb.find_one({'crbmodel.name':pluginname})
        if model == None:
            print("Invalid model name ",pluginname)
            exit
        pprint.pprint(model)
        dbid=model['uuid']
        
     
    plugin=pyfeelpp.crb.factoryCRBPlugin(pluginname)
    print("Plugin name:",plugin.name(), " dbid:",dbid);
    plugin.loadDBFromId(id=dbid,load=load);
    return plugin

def online( rbmodel, Nsampling=1 ):
    Dmu=rbmodel.parameterSpace()
    s=Dmu.sampling()
    s.sampling(Nsampling,"random")
    if ( rbmodel.isFiniteElementModelDBLoaded() ):
        rbmodel.initExporter()
    with Timer('rbmodel.run'):
        r=rbmodel.run(s.getVector())
    for x in r:
        print("mu=",x.parameter(),"s=",x.output(), " error=",x.errorBound())
        if ( rbmodel.isFiniteElementModelDBLoaded() ):            
            rbmodel.exportField("sol-"+str(x.parameter()),x)
    if ( rbmodel.isFiniteElementModelDBLoaded() ):
        rbmodel.saveExporter()      
    return r


def main(args):
     from optparse import OptionParser
     parser = OptionParser()
     parser.add_option('-N', '--Num', dest='N', help='number of parameters to evaluate the plugin', type=int, default=1)
     parser.add_option('-d', '--dir', dest='dir', help='directory location of the plugins', default="/usr/local/lib")
     parser.add_option('-n', '--name', dest='name', help='name of the plugin',type="string")
     parser.add_option('-i', '--id', dest='dbid', help='DB id to be used', type="string")
     parser.add_option('-l', '--load', dest='load', help='type of data to be loadedoaded', type="string",default="rb")
     (options, args) = parser.parse_args()
     print("members:",pyfeelpp.crb.CRBLoad.__members__)
     print("rb members:",pyfeelpp.crb.CRBLoad.__members__["rb"])

     e=pyfeelpp.core.Environment(sys.argv,pyfeelpp.crb.makeCRBOptions())
    
     model=loadPlugin(pluginname=options.name,plugindir=options.dir,load=pyfeelpp.crb.CRBLoad.__members__[options.load])
     r=online(model,options.N)
     
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
    
