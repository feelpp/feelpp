import feelpp
import feelpp.mor as mor
from feelpp.timing import *
import sys,os,time

class Online:
    """Class to execute the online phase of the reduced basis method
    """    
    def __init__(self,plugin,root):
        """construct a new Online code

        Args:
            plugin (_type_): plugin name
            root (_type_): root directory where the 
        """        
        self.crbdb = mor.CRBModelDB(pluginname=plugin,root=root)
        self.meta = crbdb.loadDBMetaData("last_modified", "")
        print(f"model: {meta.model_name} json: {meta.json_path} plugin name:{meta.plugin_name}")
        self.rbmodel=crbdb.loadDBPlugin(meta,load="rb")
        if ( self.rbmodel.isFiniteElementModelDBLoaded() ):
            rbmodel.initExporter()

    def sampling( Nsamples=1, type="random" ):
        """generate a sampling of Nsamples points in the parameter space

        Args:
            Nsamples (int, optional): number of samples. Defaults to 1.

        Returns:
            sampling: Parameter set sampling
        """        
        Dmu=self.rbmodel.parameterSpace()
        s=Dmu.sampling()
        s.sampling(Nsamples,type)
        return s

    def run( samples, export=False ):
        """
        run the model with the given samples

        Args:
            samples (sampling): sampling of the parameter set
            export (bool, optional): if True, export the reduced function in the finite element basis. Defaults to False.  

        Returns:
            np.array: the results of the run
        """        
        with Timer('rbmodel.run'):
            r=self.rbmodel.run(samples.getVector())
        for x in r:
            print("mu=",x.parameter(),"s=",x.output(), " error=",x.errorBound())
            if ( export and rbmodel.isFiniteElementModelDBLoaded() ):            
                self.rbmodel.exportField("sol-"+str(x.parameter()),x)
        if ( export and rbmodel.isFiniteElementModelDBLoaded() ):
            self.rbmodel.saveExporter()      
        return r


def main(args):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-N', '--Num', dest='N', help='number of parameters to evaluate the plugin', type=int, default=1)
    parser.add_option('-d', '--dir', dest='dir', help='directory location of the plugins', default=feelpp.Info.prefix() + '/' +feelpp.Info.libdir())
    parser.add_option('-n', '--name', dest='name', help='name of the plugin',type="string")
    parser.add_option('-i', '--id', dest='dbid', help='DB id to be used', type="string")
    parser.add_option('-l', '--load', dest='load', help='type of data to be loadedoaded', type="string",default="rb")
    (options, args) = parser.parse_args()
    print("members:",mor.CRBLoad.__members__)
    print("rb members:",mor.CRBLoad.__members__["rb"])
    print("libdir: {}".format(feelpp.Info.prefix() + '/' +feelpp.Info.libdir()))

    app=feelpp.Environment(sys.argv,mor.makeCRBOptions())
    o = Online(options.name,options.dir)
    s=o.sampling(options.N)
    res = o.run(samples=s,export=True)
    import pandas as pd
    df = pd.DataFrame(res)
    print(df)
    df.to_csv("results.csv")
     
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
    
