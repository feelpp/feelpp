# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
from os.path import dirname, basename, isfile, join
import glob

# %%

def plot_error(dirs, names, Ns, norm='l2'):
    """Plot errors obtained from many directories

    Args:
        dirs (list of str): list of path to export directories
        names (list of str): list of names for the legend
        Ns (list of int): sizes of the basis used for the computation of the errors
        norm (str, optional): Norm to plot (l2 or h1). Defaults to l2.
    """
    D = len(dirs)
    assert(D == len(names))
    NN = len(Ns)

    cs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

    fig, ax = plt.subplots(1, 1, figsize=(20,10))

    errIntL2 = 0
    errIntLinf = 0

    for d in range(D):
        mins  = np.zeros(NN)
        means = np.zeros(NN)
        maxs  = np.zeros(NN)


        for i, N in enumerate(Ns):
            df = pd.read_csv(f"{dirs[d]}/errorParams/errors{N}.csv")
            mins[i]  = df[f"{norm}(uh-uHn)"].min()
            means[i] = df[f"{norm}(uh-uHn)"].mean()
            maxs[i]  = df[f"{norm}(uh-uHn)"].max()



            errIntL2 += df[f"{norm}(uh-uH)"].mean()

            print(f"[NIRB online] {names[d]}: N={N},\n\tmin L2 error = {mins[i]},\n\tmean {norm} error = {means[i]},\n\tmax {norm} error = {maxs[i]}")

        ax.plot(Ns, mins, c=cs[d], label=f"{norm} {names[d]}")
        ax.plot(Ns, means, '--', c=cs[d])
        ax.plot(Ns, maxs, c=cs[d])
        ax.fill_between(Ns, mins, maxs, alpha=0.2, color=cs[d])


    ax.axhline(errIntL2 / D, c='k', linestyle='--', label=r"$\Vert u_H - u_h\Vert$")

    # Plot
    ax.set_title("Error " + [r"$H^1$", r"$L^2$"][norm=="l2"])
    ax.set_xlabel("N")
    ax.set_ylabel("Error")
    # ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(Ns[0], Ns[-1])
    ax.legend()



    # Save
    tikzplotlib.save("plot.tex")
    plt.show()

# plot_error('/data/scratch/saigre/feel-mbda/nirb/heat/np_1', [2, 3, 4, 5, 10, 15, 20, 25, 50, 100, 175, 200])
#plot_error('/data/home/elarif/feelppdb/nirb/heat/np_1', [1, 2, 4, 6, 10, 12, 14, 16, 20, 25, 30, 35, 40, 45, 50, 70, 80, 100])
plot_error(['/data/scratch/saigre/feel-nirb/nirb/heat/np_1/RESULTS/Rect/Greedy/', '/data/scratch/saigre/feel-nirb/nirb/heat/np_1/RESULTS/Rect/noGreedy/'],
           ['w/ Greedy', 'w/o Greedy'],
           [3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 75, 80])

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
