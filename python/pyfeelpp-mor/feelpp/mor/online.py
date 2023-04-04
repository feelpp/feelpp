import feelpp
import feelpp.mor as mor
from feelpp.timing import *
import sys, os, time
import pandas as pd

class Online:
    """Class to execute the online phase of the reduced basis method
    """
    def __init__(self, plugin, root):
        """construct a new Online code

        Args:
            plugin (str): plugin name
            root (str): root directory where the database is located
        """
        self.crbdb = mor.CRBModelDB(name=plugin, root=root)
        self.meta = self.crbdb.loadDBMetaData("last_modified", "")
        print(f"model: {self.meta.model_name}\njson: {self.meta.json_path.string()}\nplugin name:{self.meta.plugin_name}")
        self.rbmodel = self.crbdb.loadDBPlugin(self.meta, load="rb")
        if ( self.rbmodel.isFiniteElementModelDBLoaded() ):
            self.rbmodel.initExporter()

    def sampling( self, Nsamples=1, type="random" ):
        """generate a sampling of Nsamples points in the parameter space

        Args:
            Nsamples (int, optional): number of samples. Defaults to 1.

        Returns:
            sampling: Parameter set sampling
        """
        Dmu = self.rbmodel.parameterSpace()
        s = Dmu.sampling()
        s.sampling(Nsamples,type)
        return s

    def run( self, samples, export=False ):
        """
        run the model with the given samples

        Args:
            samples (sampling): sampling of the parameter set
            export (bool, optional): if True, export the reduced function in the finite element basis. Defaults to False.  

        Returns:
            np.array: the results of the run
        """
        with Timer('rbmodel.run'):
            r = self.rbmodel.run(samples.getVector())
        for x in r:
            print("mu =", x.parameter(), "s =", x.output(), "error =", x.errorBound())
            if ( export and self.rbmodel.isFiniteElementModelDBLoaded() ):            
                self.rbmodel.exportField("sol-" + str(x.parameter()), x)
        if ( export and self.rbmodel.isFiniteElementModelDBLoaded() ):
            self.rbmodel.saveExporter()
        return r


def convertToDataframe( res ):
    """Convert the result of the online phase to dataframe

    Parameters
    ----------
    res : list of feelpp.mor._mor.CRBResults
        list of results

    Returns
    -------
    pandas.DataFrame
        dataframe with the results : parameter, output, errorBound
    """
    names = res[0].parameter().parameterNames()
    df = pd.DataFrame(columns = tuple(names) + ('output', 'errorBound'))
    # df = pd.DataFrame(columns = ('parameter', 'output', 'errorBound'))
    for i, r in enumerate(res):
        p = r.parameter()
        df.loc[i] = [p.parameterNamed(p.parameterName(i)) for i in range(p.size())]+  [r.output(), r.errorBound()]
        # df.loc[i] = [r.parameter(), r.output(), r.errorBound()]
    return df


def main(args):
    libdir = os.path.join(feelpp.Info.prefix().string(), feelpp.Info.libdir().string())
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-N', '--Num', dest='N', help='number of parameters to evaluate the plugin', type=int, default=1)
    parser.add_option('-d', '--dir', dest='dir', help='directory location of the plugins', default=libdir)
    parser.add_option('-n', '--name', dest='name', help='name of the plugin',type="string")
    parser.add_option('-i', '--id', dest='dbid', help='DB id to be used', type="string")
    parser.add_option('-l', '--load', dest='load', help='type of data to be loadedoaded', type="string",default="rb")
    (options, args) = parser.parse_args()
    print("members:", mor.CRBLoad.__members__)
    print("rb members:", mor.CRBLoad.__members__["rb"])
    print("libdir:", libdir)

    app = feelpp.Environment(sys.argv, mor.makeCRBOptions())
    o = Online(options.name, options.dir)
    s = o.sampling(options.N)
    res = o.run(samples=s, export=True)
    df = convertToDataframe(res)
    print(df)
    df.to_csv("results.csv")

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
