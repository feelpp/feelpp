# %%
import pandas as pd 
# import tikzplotlib # Pour exporter le graphique
import matplotlib.pyplot as plt # Pour tracer des graphiques
import numpy as np # Pour le calcul numérique

from sklearn.linear_model import LinearRegression 
# le module scikit
from os.path import dirname, basename, isfile, join
import glob
# modules = glob.glob(join(dirname(__file__), "*.py"))
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]

# %%

def manage_dataFrame(dataDir, saveGlob=False, fileTosave=None):
    """Manage pandas data Frame (df) : conctatenate all df, get different df of the (mean, min, max)
        of all df in respect to number of snap (the number of snap is integer N such that the file
        contained in dataDir are named 'errorsN.csv')

    Args:
        dataDir (str): directory contaning all data Frame 
        saveGlob (bool): To save Global data frame in files. Defaults to False. 
        fileTosave (str) : file to save the global data frame. Defaults to None. 
    """

    files = glob.glob(join(dirname(dataDir), "*.csv"))
    NsnapList = [int(basename(f)[6:][:-4]) for f in files if isfile(f)]
    
    Ns = len(NsnapList)

    dfList = []

    for i in range(Ns):
        nk = NsnapList[i]
        file = dataDir+f"errors{nk}.csv"
        df = pd.read_csv(file, sep=',')
        df['N'] = nk 
        dfList.append(df)

    dfGlob = pd.concat(dfList)

    assert np.sum([df.shape[0] for df in dfList]) == dfGlob.shape[0]

    dfGlob.rename(columns={'l2u-uH':"l2(u-uH)", 'l2u-uHn':"l2(u-uHn)", 'lINFu-uH':"inf(u-uH)",
                         'lINFu-uHn':'inf(u-uHn)'}, inplace=True)

    
    if saveGlob and fileTosave != None :
        ## Save the global data frame in csv files
        dfGlob.to_csv(fileTosave, sep=',')
    
    lkeys = [k for k in dfGlob.keys() if k not in ['N', 'parameter']]

    # Get mean, min and max
    dfmean = dfGlob.pivot_table(values=lkeys, index='N', aggfunc=np.mean)
    dfmin = dfGlob.pivot_table(values=lkeys, index='N', aggfunc=np.min)
    dfmax = dfGlob.pivot_table(values=lkeys, index='N', aggfunc=np.max)
    
    l2 = 'l2(u-uHn)'
    inf = 'inf(u-uHn)'
    l2df = pd.DataFrame()
    l2df['Min'] = dfmin[l2]
    l2df['Max'] = dfmax[l2]
    l2df['Mean'] = dfmean[l2]
    l2df['(u-uH)'] = dfmean['l2(u-uH)']

    ## L infty norm 
    infdf = pd.DataFrame()
    infdf['Min'] = dfmin[inf]
    infdf['Max'] = dfmax[inf]
    infdf['Mean'] = dfmean[inf]
    infdf['(u-uH)'] = dfmean['inf(u-uH)']

    return l2df, infdf

# %%
def plot_dataFrame(df):
    """ plot the content of a data frame  

    Args:
        df (pandas.dataFrame): the dataFrame 
    """
    x = df.index

    for i in range(df.shape[1]):
        plt.scatter(x,df[df.keys()[i]], label=df.keys()[i])
        plt.plot(x,df[df.keys()[i]])
    
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(x[0], x[-1])
    plt.show()

def compare_dataFrame(df,dg):
    """compare the content of 2 data frame df and dg 

    Args:
        df (pandas.dataFrame): _description_
        dg (pandas.dataFrame): _description_
    """

    x = df.index
    col = ['red', 'blue', 'green']
    for i in range(df.shape[1]-1):
        plt.scatter(x, df[df.keys()[i]], c='red')
        plt.scatter(x, dg[dg.keys()[i]], c='blue')
        plt.plot(x, df[df.keys()[i]], c='red')
        plt.plot(x, dg[dg.keys()[i]], c='blue')

    plt.plot(x, df[df.keys()[-2]], c='red', label='w/O rectif')
    plt.plot(x, dg[dg.keys()[-2]], c='blue', label='w/ rectif')
    plt.scatter(x, df[df.keys()[-1]], c='green')
    plt.plot(x, dg[dg.keys()[-1]], c='green', label=r"$\Vert u_h - u_{H}\Vert$")

    
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(x[0], x[-1])
    plt.show()

# %%
if __name__ == '__main__':
    dir ="/Users/elarif2/elarif/data/feelppdb/nirb/heat/np_1/errorParams/"
    dirR ="/Users/elarif2/elarif/data/feelppdb/nirb/heat/np_1/errorParamsRectif/"
    norm = '2'
    
l2df, infdf = manage_dataFrame(dir)
l2dfR, infdfR = manage_dataFrame(dirR)


compare_dataFrame(l2df, l2dfR)
# %%
