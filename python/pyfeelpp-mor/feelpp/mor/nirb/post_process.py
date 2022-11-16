# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
from os.path import dirname, basename, isfile, join
import glob

# %% 
    # Functions to vizualise dataFrame 

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


def plot_time(dataFrame=None, csv_file=None):
    
    # assert (dataFrame==None)and(csv_file==None) 

    if csv_file==None:
        df = dataFrame
    else :
        df = pd.read_csv(csv_file)
    

    plt.plot(df['N'], df['toolbox'], '+', c='k', linestyle='--', label='Fine toolbox')
    plt.plot(df['N'], df['nirb_online'], '+-', label='NIRB w/o rect.')
    plt.plot(df['N'], df['nirb_online_rect'], '+-', label='NIRB w/ rect.')

    plt.xlabel('Number of basis function (N)')
    plt.ylabel('Time (s)')
    # plt.xscale('log')
    # plt.yscale('log')

    plt.legend()
    plt.show()
    # tikzplotlib.save("plot.tex")

def compare_time_parallel(csv_seq, csv_para):
    """compare elapsed time between parallel computing and sequential one

    Args:
        csv_seq (str): csv file of sequential datas
        csv_para (str): csv file of parallel datas 
    """
    dfseq = pd.read_csv(csv_seq)
    dfpar = pd.read_csv(csv_para)

    assert 'N' in dfseq.keys() 

    lkeys = [k for k in dfseq.keys() if k != 'N']

    s = 0 
    for k in lkeys :
        plt.figure(s)
        plt.plot(dfseq['N'], dfseq[k], label= k + ' 1 proc')
        plt.plot(dfpar['N'], dfpar[k], label= k + ' 4 proc')
        plt.legend()
        plt.grid('on')
        plt.xlabel('Number of Basis functions (N)')
        plt.ylabel('Time (s)')
        s+=1
        
    plt.show()
    # tikzplotlib.save("plot.tex")

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

def getRelativeErrors(df,h1=True):
    """Compute the relative error on a given data Frame 

    Args:
        df (pandas.dataFrame): dataFrame 
        h1 (bool, optional): if H1 norm is computed. Defaults to True.
    """

    l2keys = ['l2(uh-uHn)', 'l2(uh-uHn)rec', 'l2(uh-uhn)', 'l2(uh-uH)']
    h1keys = ['h1(uh-uHn)', 'h1(uh-uHn)rec', 'h1(uh-uhn)', 'h1(uh-uH)']
    keys = l2keys.copy() 
    if h1: keys = l2keys + h1keys 

    # df.to_numpy()

    dfRel = df[keys].copy() 
    dfRel['N']=df['N']

    for key in l2keys:
        dfRel[key] = df[key]/df['l2(uh)']
    if h1:
        for key in h1keys:
            dfRel[key]= df[key]/df['h1(uh)']

    return dfRel 


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


### Function to vizualise data Frame 
def plot_dataFrame(df, norm='L2'):
    """ plot the content of a data frame  

    Args:
        df (pandas.dataFrame): the dataFrame 
    """
    x = df.index
    keys = [k for k in df.keys() if k not in ['N', 'parameter']]

    if norm=='h1':
        nm = f"$H_1$"
    else :
        nm = f"$L_2$"

    for i in range(df.shape[1]-1):
        plt.scatter(x,df[keys[i]], label=str(keys[i]))
        plt.plot(x,df[keys[i]])
    
    plt.legend() 
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{nm} norm of Errors (in log)")
    plt.show()

def compare_dataFrams(listdf, keys='Mean', norm='l2'):

    labels =["$\lambda = 1.e^{-1}$", "$\lambda = 1.e^{-3}$", "$\lambda = 1.e^{-6}$", "$\lambda = 1.e^{-10}$", "$\lambda = 0$"]
    key_list = {'Mean':['Mean', 'Mean_rec', 'Mean_uh'], 'Max':['Max', 'Max_rec','Max_uh'],
             'Min':['Min', 'Min_rec', 'Min_uh'] }

    assert len(labels)>=len(listdf)

    if norm=='h1':
        nm = f"$H_1$"
    else :
        nm = f"$L_2$"

    idk = 1
    i=0
    for df in listdf:
        plt.scatter(df.index, df[key_list[keys][idk]], label=labels[i])
        plt.plot(df.index, df[key_list[keys][idk]])
        i+=1

    plt.scatter(df.index, df[key_list[keys][0]], label="w/o rectif")
    plt.plot(df.index, df[key_list[keys][0]])

    plt.legend()
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{nm} norm of Errors (in log scale)")
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

    normUHn = {'l2':r"$\Vert u^\mathcal{N}_h - u^N_{Hh}\Vert_{L_2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^N_{Hh}\Vert_{H_1}$"}
    normUhn = {'l2':r"$\Vert u^\mathcal{N}_h - u^N_{h}\Vert_{L_2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^N_{h}\Vert_{H_1}$"}
    normUh = {'l2':r"$\Vert u^\mathcal{N}_h - u^\mathcal{N}_{Hh}\Vert_{L_2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^\mathcal{N}_{Hh}\Vert_{H_1}$"}
  
    keyUh = {'l2': 'l2(uh-uH)', 'h1':'h1(uh-uH)'}

    if norm=='h1':
        nm = f"$H_1$"
    else :
        nm = f"$L_2$"

    plt.scatter(xf, df[key_list[keys][0]], c='red', label=normUHn[norm] + 'w/o rectif')
    plt.scatter(xf, df[key_list[keys][1]], c='blue', label=normUHn[norm] + 'w/ rectif')
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
    plt.ylabel(f"{nm} norm of Errors (in log scale)")
    # tikzplotlib.save(f"{keys}.tex")
    plt.show()


def compare_2dataStats(df,dg, keys='Mean', norm='l2'):
    """compare two dataFrames in respect of given keys and norm 
        keys take 'Min', 'Max' or 'Mean' 

    Args:
        df (pandas.dataFrame): _description_
        dg (pandas.dataFrame): _description_
        keys (str) : Min, Max or Mean 
    """

    xf = df.index
    xg = dg.index 

    key_list = {'Mean':['Mean', 'Mean_rec', 'Mean_uh'], 'Max':['Max', 'Max_rec','Max_uh'],
             'Min':['Min', 'Min_rec', 'Min_uh'] }

    normUHn = {'l2':r"$\Vert u^\mathcal{N}_h - u^N_{Hh}\Vert_{L_2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^N_{Hh}\Vert_{H_1}$"}
    normUhn = {'l2':r"$\Vert u^\mathcal{N}_h - u^N_{h}\Vert_{L_2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^N_{h}\Vert_{H_1}$"}
    normUh = {'l2':r"$\Vert u^\mathcal{N}_h - u^\mathcal{N}_{Hh}\Vert_{L_2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^\mathcal{N}_{Hh}\Vert_{H_1}$"}

    keyUh = {'l2': 'l2(uh-uH)', 'h1':'h1(uh-uH)'}

    if norm=='h1':
        nm = f"$H_1$"
    else :
        nm = f"$L_2$"

    plt.scatter(xf, df[key_list[keys][0]], marker='o', c='red', label=normUHn[norm] + ' w/o rect - $\mathbb{P}_1$')
    plt.scatter(xg, dg[key_list[keys][0]], marker='x', c='red', label=normUHn[norm] + ' w/o rect - $\mathbb{P}_2$')

    plt.scatter(xf, df[key_list[keys][1]], marker='o', c='blue', label=normUHn[norm] + ' w/ rectif - $\mathbb{P}_1$')
    plt.scatter(xg, dg[key_list[keys][1]], marker='x', c='blue', label=normUHn[norm] + ' w/ rectif - $\mathbb{P}_2$')

    plt.plot(xf, df[key_list[keys][0]], c='red')
    plt.plot(xg, dg[key_list[keys][0]], c='red')

    plt.plot(xf, df[key_list[keys][1]], c='blue')
    plt.plot(xg, dg[key_list[keys][1]], c='blue')
    
    plt.legend()
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{nm} norm of Errors (in log scale)")
    # tikzplotlib.save(f"{keys}.tex")
    plt.show()




# %%
### Main 

if __name__ == "__main__":
    
    import sys 

    # Warning to specifie this env path !!!!
    envpath = "/Users/elarif2/elarif/devel/docker.feel/feelppdb/nirb/heat/np_1/" 

    #%% Get dataFrame from csv file 
    file = envpath + "errors50ParamsLambda10P1.csv"
    # load absolute errors
    dfGlob = pd.read_csv(file, sep=',')
    # compute relative errors 
    dfRel = getRelativeErrors(dfGlob)

    #%%
    ### Get stats for all parameters 

    l2df, h1df   = getDataStat(dfGlob) # l1 and h1 error associated 
    l2dfRel, h1dfRel = getDataStat(dfRel) # // 

    plot_dataFrame(l2df)
    compare_dataStats(l2df, keys='Mean')

    # %%
    ## Tronctae error datas into N parameter (N<=50)
    N = 5
    dfN = troncateNparam(dfGlob, N) 
    
    l2dfN, h1dfN   = getDataStat(dfN) # l1 and h1 error associated with statistical infos 
    
    plot_dataFrame(l2dfN)
    compare_dataStats(l2dfN, keys='Mean')

    #%%
    # Compare rectification according to regularization parameter (\lambda) 
    shortfiles = ["errors50ParamsLambda1P1.csv", "errors50ParamsLambda3P1.csv","errors50ParamsLambda6P1.csv", 
                    "errors50ParamsLambda10P1.csv", "errors50ParamsLambda0P1.csv"]
    
    listdfl2, listdfh1 = [], []

    for st in shortfiles:
        file = envpath + st 
        dfG = pd.read_csv(file, sep=',')
        dl2, dh1 = getDataStat(dfG)
        listdfl2.append(dl2)
        listdfh1.append(dh1)
    
    compare_dataFrams(listdfl2)

    # %%
    # Compare P1 vs P2 solvers  
    filep1 = envpath + "errors50ParamsLambda10P1.csv"
    filep2 = envpath + "errors50ParamsLambda10P2.csv"

    dfp1 = pd.read_csv(filep1, sep=',')
    dfRelp1 = getRelativeErrors(dfp1)

    dfp2 = pd.read_csv(filep2, sep=',')
    dfRelp2 = getRelativeErrors(dfp2)

    l2p1, h1p1 = getDataStat(dfp1)
    l2p2, h1p2 = getDataStat(dfp2)

    compare_2dataStats(l2p1, l2p2)


    # %%
    ## Compare execution time between toolbox and nirb online w/ and w/o rectification  
    file = envpath + "nirbOnline_time_exec_norect.dat"
    npn = np.loadtxt(file)
    file =  envpath + "nirbOnline_time_exec_rect.dat"
    npr = np.loadtxt(file)
    dic = {}
    dic['N'] = npn[:,0]
    dic['toolbox'] = npr[:,1]
    dic['nirb_online'] = npn[:,2]
    dic['nirb_online_rect']= npr[:,2]
    df = pd.DataFrame(dic)
    
    plot_time(df)

    #%% 
    # Compare execution time  between parallel and sequential  
    file1 = envpath + "nirb_time_exec_rect.csv"
    file4 = envpath[:-5] + "np_4/nirb_time_exec_rect.csv"

    compare_time_parallel(file1, file4) 

