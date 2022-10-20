# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
from os.path import dirname, basename, isfile, join
import glob

# %%

def plot_error(dir_name, Ns):
    N = len(Ns)
    # l2u-uHn, lINFu-uHn, l2u-uH, lINFu-uH
    mins2Nirb_norect = np.zeros(N)
    maxs2Nirb_norect = np.zeros(N)
    moy2Nirb_norect  = np.zeros(N)

    minsInfNirb_norect = np.zeros(N)
    maxsInfNirb_norect = np.zeros(N)
    moyInfNirb_norect  = np.zeros(N)

    mins2Int_norect = np.zeros(N)
    maxs2Int_norect = np.zeros(N)
    moy2Int_norect  = np.zeros(N)

    minsInfInt_norect = np.zeros(N)
    maxsInfInt_norect = np.zeros(N)
    moyInfInt_norect  = np.zeros(N)


    mins2Nirb_rect = np.zeros(N)
    maxs2Nirb_rect = np.zeros(N)
    moy2Nirb_rect  = np.zeros(N)

    minsInfNirb_rect = np.zeros(N)
    maxsInfNirb_rect = np.zeros(N)
    moyInfNirb_rect  = np.zeros(N)

    mins2Int_rect = np.zeros(N)
    maxs2Int_rect = np.zeros(N)
    moy2Int_rect  = np.zeros(N)

    minsInfInt_rect = np.zeros(N)
    maxsInfInt_rect = np.zeros(N)
    moyInfInt_rect  = np.zeros(N)

    for i,N in enumerate(Ns):
        df_norect = pd.read_csv(f"{dir_name}/errorParams/errors{N}.csv")
        mins2Nirb_norect[i] = df_norect["l2u-uHn"].min()
        maxs2Nirb_norect[i] = df_norect["l2u-uHn"].max()
        moy2Nirb_norect[i]  = df_norect["l2u-uHn"].mean()

        minsInfNirb_norect[i] = df_norect["lINFu-uHn"].min()
        maxsInfNirb_norect[i] = df_norect["lINFu-uHn"].max()
        moyInfNirb_norect[i]  = df_norect["lINFu-uHn"].mean()

        mins2Int_norect[i] = df_norect["l2u-uH"].min()
        maxs2Int_norect[i] = df_norect["l2u-uH"].max()
        moy2Int_norect[i]  = df_norect["l2u-uH"].mean()

        minsInfInt_norect[i] = df_norect["lINFu-uH"].min()
        maxsInfInt_norect[i] = df_norect["lINFu-uH"].max()
        moyInfInt_norect[i]  = df_norect["lINFu-uH"].mean()

        df_rect = pd.read_csv(f"{dir_name}/errorParamsRectif/errors{N}.csv")
        mins2Nirb_rect[i] = df_rect["l2u-uHn"].min()
        maxs2Nirb_rect[i] = df_rect["l2u-uHn"].max()
        moy2Nirb_rect[i]  = df_rect["l2u-uHn"].mean()

        minsInfNirb_rect[i] = df_rect["lINFu-uHn"].min()
        maxsInfNirb_rect[i] = df_rect["lINFu-uHn"].max()
        moyInfNirb_rect[i]  = df_rect["lINFu-uHn"].mean()

        mins2Int_rect[i] = df_rect["l2u-uH"].min()
        maxs2Int_rect[i] = df_rect["l2u-uH"].max()
        moy2Int_rect[i]  = df_rect["l2u-uH"].mean()

        minsInfInt_rect[i] = df_rect["lINFu-uH"].min()
        maxsInfInt_rect[i] = df_rect["lINFu-uH"].max()
        moyInfInt_rect[i]  = df_rect["lINFu-uH"].mean()

    print("MinNitb_norect", mins2Nirb_norect)
    print("MaxNitb_norect", maxs2Nirb_norect)
    print("MoyNitb_norect", moy2Nirb_norect)

    print("MinInfNitb_norect", minsInfNirb_norect)
    print("MaxInfNitb_norect", maxsInfNirb_norect)
    print("MoyInfNitb_norect", moyInfNirb_norect)

    print("MinInt_rect", mins2Int_rect)
    print("MaxInt_rect", maxs2Int_rect)
    print("MoyInt_rect", moy2Int_rect)

    print("MinInfInt_rect", minsInfInt_rect)
    print("MaxInfInt_rect", maxsInfInt_rect)
    print("MoyInfInt_rect", moyInfInt_rect)


    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(20,10))

    ax[0].plot(Ns, mins2Nirb_norect, c="tab:blue", label=r"$\Vert u_h - u_{hH}^N\Vert_{L^2}$ w/o rect.")
    ax[0].plot(Ns, maxs2Nirb_norect, c="tab:blue")
    ax[0].plot(Ns, moy2Nirb_norect, c="tab:blue")
    ax[0].fill_between(Ns, mins2Nirb_norect, maxs2Nirb_norect, alpha=0.2)
    ax[0].plot(Ns, mins2Nirb_rect, c="tab:orange", label=r"$\Vert u_h - u_{hH}^N\Vert_{L^2}$ w/ rect.")
    ax[0].plot(Ns, maxs2Nirb_rect, c="tab:orange")
    ax[0].plot(Ns, moy2Nirb_rect, c="tab:orange")
    ax[0].fill_between(Ns, mins2Nirb_rect, maxs2Nirb_rect, alpha=0.2)
    ax[0].axhline(mins2Int_norect[0], color="tab:green", linestyle="--", label=r"$\Vert u_h - u_H\Vert_{L^2}$")
    ax[0].axhline(maxs2Int_norect[0], color="tab:green", linestyle="--")
    ax[0].axhline(moy2Int_norect[0], color="tab:green", linestyle="--")
    ax[0].set_title(r"Error $L^2$")
    ax[0].set_xlabel("N")
    ax[0].set_ylabel("Error")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(Ns[0], Ns[-1])
    ax[0].legend()

    ax[1].plot(Ns, minsInfNirb_norect, c="tab:blue", label=r"$\Vert u_h - u_{hH}^N\Vert_{\infty}$ w/o rect.")
    ax[1].plot(Ns, maxsInfNirb_norect, c="tab:blue")
    ax[1].plot(Ns, moyInfNirb_norect, c="tab:blue")
    ax[1].fill_between(Ns, minsInfNirb_norect, maxsInfNirb_norect, alpha=0.2)
    ax[1].plot(Ns, minsInfNirb_rect, c="tab:orange", label=r"$\Vert u_h - u_{hH}^N\Vert_{\infty}$ w/ rect.")
    ax[1].plot(Ns, maxsInfNirb_rect, c="tab:orange")
    ax[1].plot(Ns, moyInfNirb_rect, c="tab:orange")
    ax[1].fill_between(Ns, minsInfNirb_rect, maxsInfNirb_rect, alpha=0.2)
    ax[1].axhline(minsInfInt_norect[0], color="tab:green", linestyle="--", label=r"$\Vert u_h - u_{hH}^N\Vert_{\infty}$")
    ax[1].axhline(maxsInfInt_norect[0], color="tab:green", linestyle="--")
    ax[1].axhline(moyInfInt_norect[0], color="tab:green", linestyle="--")
    ax[1].set_title(r"Error $L^\infty$")
    ax[1].set_xlabel("N")
    ax[1].set_ylabel("Error")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(Ns[0], Ns[-1])
    ax[1].legend()

    # Save
    tikzplotlib.save("plot.tex")
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


    # dir ="/Users/elarif2/elarif/data/feelppdb/nirb/heat/np_1/errorParams50/"
    # dirR ="/Users/elarif2/elarif/data/feelppdb/nirb/heat/np_1/errorParamsRectif50/"

    dir  = str(sys.argv[1]) # directory containing errors w/O rectification 
    dirR = str(sys.argv[2]) # // w/ rectification 
    norm ='l2'

    dfGlob  = concat_dataFrame(dir) # get global dataFrame w/O rectification 
    dfGlobR = concat_dataFrame(dirR) # global dataFrame w/ rectification 

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
    l2df, h1df   = getDataStat(dkN) # l1 and h1 error associated 
    l2dfR, h1dfR = getDataStat(dkNR) # // 

    plot_dataFrame(l2df)
    compare_dataFrame(l2df, l2dfR, keys='Max')