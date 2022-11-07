# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import tikzplotlib
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
def getDataStat(df, h1=True):
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
    
    l2 = 'l2(uh-uHn)'
    l2rec = 'l2(uh-uHn)rec'
    l2uh = 'l2(uh-uhn)'

    h1rec = 'h1(uh-uHn)rec'
    h1 = 'h1(uh-uHn)'
    h1uh = 'h1(uh-uhn)'

    
    ## L2 norm 
    l2df = pd.DataFrame()

    l2df['Min'] = dfmin[l2]
    l2df['Min_rec'] = dfmin[l2rec]
    l2df['Min_uh'] = dfmin[l2uh]

    l2df['Max'] = dfmax[l2]
    l2df['Max_rec'] = dfmax[l2rec]
    l2df['Max_uh'] = dfmax[l2uh]
    
    l2df['Mean'] = dfmean[l2]
    l2df['Mean_rec'] = dfmean[l2rec]
    l2df['Mean_uh'] = dfmean[l2uh]

    l2df['l2(uh-uH)'] = dfmean['l2(uh-uH)']

    ## H1 norm
    h1df = pd.DataFrame() 
    if h1:

        h1df['Min'] = dfmin[h1]
        h1df['Min_rec'] = dfmin[h1rec]
        h1df['Min_uh'] = dfmin[h1uh]

        h1df['Max'] = dfmax[h1]
        h1df['Max_rec'] = dfmax[h1rec]
        h1df['Max_uh'] = dfmax[h1uh]
        
        h1df['Mean'] = dfmean[h1]
        h1df['Mean_rec'] = dfmean[h1rec]
        h1df['Mean_uh'] = dfmean[h1uh]

        h1df['h1(uh-uH)'] = dfmean['h1(uh-uH)']

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
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{norm} norm of Errors (in log)")
    plt.show()

def compare_dataStats(df, keys='Mean', norm='l2'):
    """plot statiscal comparison of a dataFrame in respect of given keys 
        keys take 'Min', 'Max' or 'Mean'

        This will plot computed errors with and without rectification 

    Args:
        df (pandas.dataFrame): _description_
        dg (pandas.dataFrame): _description_
        keys (str) : Min, Max or Mean 
    """

    xf = df.index

    key_list = {'Mean':['Mean', 'Mean_rec', 'Mean_uh'], 'Max':['Max', 'Max_rec','Max_uh'],
             'Min':['Min', 'Min_rec', 'Min_uh'] }

    normUHn = {'l2':r"$\Vert u_h - u_{NH}\Vert_{L_2}$" , 'h1': r"$\Vert u_h - u_{NH}\Vert_{H_1}$"}
    normUhn = {'l2':r"$\Vert u_h - u_{Nh}\Vert_{L_2}$" , 'h1': r"$\Vert u_h - u_{Nh}\Vert_{H_1}$"}
    normUh = {'l2':r"$\Vert u_h - u_{H}\Vert_{L_2}$" , 'h1': r"$\Vert u_h - u_{H}\Vert_{H_1}$"}
    keyUh = {'l2': 'l2(uh-uH)', 'h1':'h1(uh-uH)'}

    plt.scatter(xf, df[key_list[keys][0]], c='red', label=normUHn[norm] + 'w/o rectif')
    plt.scatter(xf, df[key_list[keys][1]], c='blue', label=normUHn[norm] + 'w/o rectif')
    plt.scatter(xf, df[key_list[keys][2]], c='green', label=normUhn[norm])

    plt.plot(xf, df[key_list[keys][0]], c='red')
    plt.plot(xf, df[key_list[keys][1]], c='blue')
    plt.plot(xf, df[key_list[keys][2]], c='green')

    # interpolate norm 
    plt.scatter(xf, df[keyUh[norm]], label=normUh[norm])
    plt.plot(xf, df[keyUh[norm]])

    
    plt.legend()
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{norm} norm of Errors (in log scale)")
    plt.show()



# %%


if __name__ == "__main__":
    import sys
    # plot_error(sys.argv[1], sys.argv[2:])
    # plot_time(sys.argv[1])

    ### Read global datas 
    # file = "/data/home/elarif/nirbDatas/heat/errors50Params.csv"  # From gaya 
    # fileR = "/data/home/elarif/nirbDatas/heat/errors50ParamsRectif.csv"
    
    # file  = str(sys.argv[1]) # file containing errors w/O rectification 
    # fileR = str(sys.argv[2]) # // w/ rectification 

    # file = "/feel/feelppdb/nirb/heat/np_1/errorRelative.csv"
    file = "/Users/elarif2/elarif/devel/docker.feel/feelppdb/nirb/heat/np_1/errors50Params.csv"
    # norm ='l2'

    dfGlob = pd.read_csv(file, sep=',')
    # dfGlobR = pd.read_csv(fileR, sep=',')


    # %%   
    ### Get stats for 50 parameters 

    l2df, h1df   = getDataStat(dfGlob) # l1 and h1 error associated 
    # l2dfR, h1dfR = getDataStat(dfGlobR) # // 

    # plot_dataFrame(l2df)
    compare_dataStats(l2df, keys='Mean')

    # %%
    
    # file = "/feel/feelppdb/nirb/heat/np_1/errorRelative.csv"
    # file = "/Users/elarif2/elarif/devel/docker.feel/feelppdb/nirb/heat/np_1/errors50Params.csv"
    # norm ='l2'

    # dfGlob = pd.read_csv(file, sep=',')
    # plot_dataFrame(dfGlob)

    # %%
    ## Tronctae into N parameter (N<=50)

    # dkN = troncateNparam(dfGlob, 1) 
    # dkNR = troncateNparam(dfGlobR, 1)

    # l2df, h1df   = getDataStat(dkN) # l1 and h1 error associated with statistical infos 
    # l2dfR, h1dfR = getDataStat(dkNR) # // 

    # plot_dataFrame(l2df)
    # compare_dataFrame(l2df, l2dfR, keys='Max', norm=norm)