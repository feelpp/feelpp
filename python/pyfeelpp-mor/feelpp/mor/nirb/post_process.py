# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
from os.path import dirname, basename, isfile, join
import glob

# %%

def plot_error(dirs, names, Ns):
    """Plot errors obtained from many directories

    Args:
        dirs (list of str): list of path to export directories
        names (list of str): list of names for the legend
        Ns (list of int): sizes of the basis used for the computation of the errors
    """
    D = len(dirs)
    assert(D == len(names))
    NN = len(Ns)

    cs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

    fig, ax = plt.subplots(1, 2, figsize=(20,10))

    errIntL2 = 0
    errIntLinf = 0

    for d in range(D):
        minsL2  = np.zeros(NN)
        meansL2 = np.zeros(NN)
        maxsL2  = np.zeros(NN)

        minsLinf  = np.zeros(NN)
        meansLinf = np.zeros(NN)
        maxsLinf  = np.zeros(NN)

        for i, N in enumerate(Ns):
            df = pd.read_csv(f"{dirs[d]}/errorParams/errors{N}.csv")
            minsL2[i]  = df["l2u-uHn"].min()
            meansL2[i] = df["l2u-uHn"].mean()
            maxsL2[i]  = df["l2u-uHn"].max()

            minsLinf[i]  = df["lINFu-uHn"].min()
            meansLinf[i] = df["lINFu-uHn"].mean()
            maxsLinf[i]  = df["lINFu-uHn"].max()

            errIntL2 += df["l2u-uH"].mean()
            errIntLinf += df["lINFu-uH"].mean()

            print(f"[NIRB online] {names[d]}: N={N},\n\tmin L2 error = {minsL2[i]},\n\tmean L2 error = {meansL2[i]},\n\tmax L2 error = {maxsL2[i]}")
            print(f"[NIRB online] {names[d]}: N={N},\n\tmin Linf error = {minsLinf[i]},\n\tmean Linf error = {meansLinf[i]},\n\tmax Linf error = {maxsLinf[i]}")

        ax[0].plot(Ns, minsL2, c=cs[d], label=f"L2 {names[d]}")
        ax[0].plot(Ns, meansL2, '--', c=cs[d])
        ax[0].plot(Ns, maxsL2, c=cs[d])
        ax[0].fill_between(Ns, minsL2, maxsL2, alpha=0.2, color=cs[d])

        ax[1].plot(Ns, minsLinf, c=cs[d], label=f"Linf {names[d]}")
        ax[1].plot(Ns, meansLinf, '--', c=cs[d])
        ax[1].plot(Ns, maxsLinf, c=cs[d])
        ax[1].fill_between(Ns, minsLinf, maxsLinf, alpha=0.2, color=cs[d])

    ax[0].axhline(errIntL2 / D, c='k', linestyle='--', label=r"$\Vert u_H - u_h\Vert$")
    ax[1].axhline(errIntLinf / D, c='k', linestyle='--', label=r"$\Vert u_H - u_h\Vert$")

    # Plot
    ax[0].set_title(r"Error $L^2$")
    ax[0].set_xlabel("N")
    ax[0].set_ylabel("Error")
    # ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(Ns[0], Ns[-1])
    ax[0].legend()

    ax[1].set_title(r"Error $L^\infty$")
    ax[1].set_xlabel("N")
    ax[1].set_ylabel("Error")
    # ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(Ns[0], Ns[-1])
    ax[1].legend()

    # Save
    # tikzplotlib.save("plot.tex")
    plt.show()

# plot_error('/data/scratch/saigre/feel-mbda/nirb/heat/np_1', [2, 3, 4, 5, 10, 15, 20, 25, 50, 100, 175, 200])
#plot_error('/data/home/elarif/feelppdb/nirb/heat/np_1', [1, 2, 4, 6, 10, 12, 14, 16, 20, 25, 30, 35, 40, 45, 50, 70, 80, 100])


# %%
def plot_time(csv_file):
    df = pd.read_csv(csv_file)
    plt.axhline(df['time_toolbox'].mean(), c='k', linestyle='--', label='Fine toolbox')
    plt.plot(df['N'], df['time_nirb'], '+-', label='NIRB w/o rect.')
    plt.plot(df['N'], df['time_nirb_rect'], '+-', label='NIRB w/ rect.')

    plt.xlabel('N')
    plt.ylabel('Time (s)')
    # plt.xscale('log')
    plt.yscale('log')

    plt.legend()
    plt.show()
    # tikzplotlib.save("plot.tex")


# %%
### Manage data Frame 

def concat_dataFrame(dataDir, saveGlob=False, fileTosave=None):
    """ conctatenate all df from directory dataDir in a one single dataFrame
        by adding a series of the number of snapshot associated (the number of snap is integer N such that the file
        contained in dataDir are named 'errorsN.csv')

    Args:
        dataDir (str): directory contaning all data Frame 
        saveGlob (bool): To save Global data frame in files. Defaults to False. 
        fileTosave (str) : file to save the global data frame. Defaults to None. 
    """

    files = glob.glob(join(dirname(dataDir), "*.csv"))
    NsnapList = [int(basename(f)[6:][:-4]) for f in files if isfile(f)]
    
    Ns = len(NsnapList)

    nk = NsnapList[0]
    file = dataDir+f"errors{nk}.csv"
    dfGlob = pd.read_csv(file, sep=',')
    dfGlob['N'] = nk 
    allSize = dfGlob.shape[0]
    
    for i in range(1,Ns):
        nk = NsnapList[i]
        file = dataDir+f"errors{nk}.csv"
        df = pd.read_csv(file, sep=',')
        df['N'] = nk
        allSize += df.shape[0]
        dfGlob = pd.concat([dfGlob, df])

    assert allSize == dfGlob.shape[0], f" size sum = {allSize}, size Glob = {dfGlob.shape[0]}"

    dfGlob.rename(columns={'l2u-uH':"l2(u-uH)", 'l2u-uHn':"l2(u-uHn)", 'lINFu-uH':"inf(u-uH)",
                         'lINFu-uHn':'inf(u-uHn)'}, inplace=True)

    
    if saveGlob and fileTosave != None :
        ## Save the global data frame in csv files
        dfGlob.to_csv(fileTosave, sep=',')
    
    return dfGlob


def getDataStat(df):
    """ Get some statistic infos from a given dataFrame df : 
        (Min, Max, Mean) in respect to number of snapshot series 

    Args:
        df (pandas.DataFrame) : the global data frame 
    
    Return :
        l2df(pandas.DataFrame) : data frame associated to stat relative to L2 norm
        h1df(pandas.DataFrame) : data frame associated to stat relative to H1 norm
    """

    lkeys = [k for k in df.keys() if k not in ['N', 'parameter']]

    # Get mean, min and max
    dfmean = df.pivot_table(values=lkeys, index='N', aggfunc=np.mean)
    dfmin = df.pivot_table(values=lkeys, index='N', aggfunc=np.min)
    dfmax = df.pivot_table(values=lkeys, index='N', aggfunc=np.max)
    
    l2 = 'l2(u-uHn)'
    h1 = 'h1(u-uHn)'
    ## L2 norm 
    l2df = pd.DataFrame()
    l2df['Min'] = dfmin[l2]
    l2df['Max'] = dfmax[l2]
    l2df['Mean'] = dfmean[l2]
    l2df['l2(u-uH)'] = dfmean['l2(u-uH)']

    ## H1 norm 
    h1df = pd.DataFrame()
    h1df['Min'] = dfmin[h1]
    h1df['Max'] = dfmax[h1]
    h1df['Mean'] = dfmean[h1]
    h1df['h1(u-uH)'] = dfmean['h1(u-uH)']
    

    return l2df, h1df

def troncateNparam(df, Nparam=1, start='first'):
    """ get errors for given number of parameter (Nparam)
            Nparam <= 50
    Args:
        df (pandas.dataFrame): data frame
        Nparam (int) : number of parameter to troncate from data frame. Defaults to 1. 
        start (srt) : starting trocate for first, or last element of sampling of 50 parameters
    """

    assert Nparam <= 50, f"Nparam should be <= 50, Nparam = {Nparam}"

    listN = list(df['N'])
    listN = list(set(listN))
    listkeys = [k for k in df.keys() if k != 'parameter']
    ind = np.where(df['N']==listN[0])[0][:Nparam]

    dg = df.iloc[ind,:]
    dg = dg[listkeys]
    dk = pd.DataFrame(dict(dg))
    for i in listN[1:]:
        ind = np.where(df['N']==i)[0][:Nparam]
        dg =  df.iloc[ind,:]
        dg = dg[listkeys]
        dk = pd.concat([dk, pd.DataFrame(dict(dg))])

    # dk.sort_values('N', inplace=True)
    dk.index = range(dk.shape[0])
    return dk 

# %%

### Vizualise data Frame 
def plot_dataFrame(df, norm='L2'):
    """ plot the content of a data frame  

    Args:
        df (pandas.dataFrame): the dataFrame 
    """
    x = df.index 
    keys = [k for k in df.keys() if k not in ['N', 'parameter']]

    for i in range(df.shape[1]-1):
        plt.scatter(x,df[keys[i]], label=str(keys[i]))
        plt.plot(x,df[keys[i]])
    
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{norm} norm of Errors (in log)")
    plt.xlim(x[0], x[-1])
    plt.show()

def compare_dataFrame(df,dg, keys='Mean', norm='l2'):
    """compare 2 dataFrame in espect of given keys 
        keys take 'Min', 'Max' or 'Mean'

    Args:
        df (pandas.dataFrame): _description_
        dg (pandas.dataFrame): _description_
        keys (str) : Min, Max or Mean 
    """

    xf = df.index
    xg = dg.index

    normUhn = {'l2':r"$\Vert u_h - u_{Hh}\Vert_{L_2}$" , 'h1': r"$\Vert u_h - u_{Hh}\Vert_{H_1}$"}
    normUh = {'l2':r"$\Vert u_h - u_{H}\Vert_{L_2}$" , 'h1': r"$\Vert u_h - u_{H}\Vert_{H_1}$"}
    keyUh = {'l2': 'l2(u-uH)', 'h1':'h1(u-uH)'}

    plt.scatter(xf, df[keys], c='red')
    plt.scatter(xg, dg[keys], c='blue')
    plt.plot(xf, df[keys], c='red')
    plt.plot(xg, dg[keys], c='blue')

    plt.plot(xf, df[keys], c='red', label=normUhn[norm] + 'w/O rectif')
    plt.plot(xg, dg[keys], c='blue', label=normUhn[norm] + 'w/ rectif')

    plt.scatter(xf, df[keyUh[norm]], c='green')
    plt.plot(xg, dg[keyUh[norm]], c='green', label=normUh[norm])

    
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{norm} norm of Errors (in log)")
    plt.xlim(xf[0], xf[-1])
    plt.show()



# %%


if __name__ == "__main__":
    import sys
    # plot_error(sys.argv[1], sys.argv[2:])
    # plot_time(sys.argv[1])

    ### Read global datas 
    file = "/data/home/elarif/nirbDatas/heat/errors50Params.csv"  # From gaya 
    fileR = "/data/home/elarif/nirbDatas/heat/errors50ParamsRectif.csv"
    
    # file  = str(sys.argv[1]) # file containing errors w/O rectification 
    # fileR = str(sys.argv[2]) # // w/ rectification 
    norm ='l2'

    dfGlob = pd.read_csv(file, sep=',')
    dfGlobR = pd.read_csv(fileR, sep=',')


    # %%   
    ### Get stats for 50 parameters 

    l2df, h1df   = getDataStat(dfGlob) # l1 and h1 error associated 
    l2dfR, h1dfR = getDataStat(dfGlobR) # // 

    plot_dataFrame(l2df)
    compare_dataFrame(l2df, l2dfR, keys='Min')

    # %%
    ## Tronctae into N parameter (N<=50)

    dkN = troncateNparam(dfGlob, 1) 
    dkNR = troncateNparam(dfGlobR, 1)

    l2df, h1df   = getDataStat(dkN) # l1 and h1 error associated with statistical infos 
    l2dfR, h1dfR = getDataStat(dkNR) # // 

    plot_dataFrame(l2df)
    compare_dataFrame(l2df, l2dfR, keys='Max', norm=norm)
