import pandas as pd
pd.options.plotting.backend = "plotly"
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
import os
import tikzplotly

## Functions to vizualise dataFrame

def plotRIC(RIC):
    """plot the RIC value in respect to basis function

    Parameters
    ----------
    RIC : list or numpy.array
        tab of ric value
    """
    plt.plot(np.arange(RIC.size), RIC)
    plt.title("RIC")
    plt.xlabel("Basis function")
    plt.ylabel("RIC values")

def plotTime(dataFrame=None, csv_file=None):
    """Plot execution time of nirb method given either in dataFrame or in csv file

    Args:
    -----
        dataFrame (pandas.dataFrame, optional): dataFrame of the execution time. Defaults to None.
        csv_file (str, optional): csv file of execution itme. Defaults to None.
    """

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
    tikzplotlib.save("plotTime.tex")

def compareTime(csv1, csv2, type='nirb_offline'):
    """Compare elapsed time between two methods

    Args:
    ----
        csv1 (str): csv file of first method
        csv2 (str): csv file of second method
        type (str): type of result to plot. Defaults to 'nirb_offline' ('nirb_online', 'toolbox')
    """
    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)

    assert 'N' in df1.keys()
    assert type in df1.keys(), f'type not in keys'


    plt.figure()
    plt.plot(df1['N'], df1[type], label= type + ' w/o Greedy')
    plt.plot(df2['N'], df2[type], label= type + ' w/ Greedy')
    plt.legend()
    plt.grid('on')
    plt.xlabel('Number of Basis functions (N)')
    plt.ylabel('Time (s)')

    plt.show()
    # tikzplotlib.save("compareTime.tex")

### Manage data Frame
def getDataStat(df, h1norm=True):
    """Get some statistic infos from a given dataFrame df :
        (Min, Max, Mean) in respect to number of snapshot series

    Args:
    -----
        df (pandas.DataFrame) : the global data frame
        h1 (bool) : to get statistical infos in H1 norm

    Returns:
    --------
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

    l2df['l2NirbMin'] = dfmin[l2]
    l2df['l2NirbMean'] = dfmean[l2]
    l2df['l2NirbMax'] = dfmax[l2]

    l2df['l2NirbRectMin'] = dfmin[l2rec]
    l2df['l2NirbRectMean'] = dfmean[l2rec]
    l2df['l2NirbRectMax'] = dfmax[l2rec]

    l2df['l2ProjMin'] = dfmin[l2uh]
    l2df['l2ProjMean'] = dfmean[l2uh]
    l2df['l2ProjMax'] = dfmax[l2uh]

    l2df['l2IntMin'] = dfmin['l2(uh-uH)']
    l2df['l2IntMean'] = dfmean['l2(uh-uH)']
    l2df['l2IntMax'] = dfmax['l2(uh-uH)']


    ## H1 norm
    if h1norm:

        h1df = pd.DataFrame()

        h1df['h1NirbMin'] = dfmin[h1]
        h1df['h1NirbMean'] = dfmean[h1]
        h1df['h1NirbMax'] = dfmax[h1]

        h1df['h1NirbRectMin'] = dfmin[h1rec]
        h1df['h1NirbRectMean'] = dfmean[h1rec]
        h1df['h1NirbRectMax'] = dfmax[h1rec]

        h1df['h1ProjMin'] = dfmin[h1uh]
        h1df['h1ProjMean'] = dfmean[h1uh]
        h1df['h1ProjMax'] = dfmax[h1uh]

        h1df['h1IntMin'] = dfmin['h1(uh-uH)']
        h1df['h1IntMean'] = dfmean['h1(uh-uH)']
        h1df['h1IntMax'] = dfmax['h1(uh-uH)']

        return l2df, h1df
    else :
        return l2df


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
        the dataFrame df should be the result of function getDataStat()
    Args:
    -----
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


### Functions to vizualise data Frame
def plotDataFrame(df, norm='l2', texSave=False):
    """ plot the content of a data frame

    Args:
    -----
        df (pandas.dataFrame): the dataFrame
        norm (str) : which norm to display
        texSave (bool) : save tex file or not. Defaults to False.
    """
    x = df.index
    keys = [k for k in df.keys() if k not in ['N', 'parameter']]

    if norm=='h1':
        nm = f"$H^1$"
    else :
        nm = f"$L^2$"

    for i in range(df.shape[1]-1):
        plt.scatter(x, df[keys[i]], label=str(keys[i]))
        plt.plot(x, df[keys[i]])

    plt.legend()
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{nm} norm of Errors (in log)")
    if texSave :
        tikzplotlib.save(f"plotError{keys}.tex")
    plt.show()

def compareListOfDataFrams(listdf, keys='Mean', norm='l2'):
    """Compare some data Frame containing in a list the dataFrames should be the result of function getDataStat()

    Args:
    -----
        listdf (list): list of data Frams
        keys (str, optional): statistic to compare (Mean, Max, Min). Defaults to 'Mean'.
        norm (str, optional): type of norm to compare. Defaults to 'l2'.
    """

    labels =["$\lambda = 1.e^{-1}$", "$\lambda = 1.e^{-3}$", "$\lambda = 1.e^{-6}$", "$\lambda = 1.e^{-10}$", "$\lambda = 0$"]
    key_list = {'Mean':['Mean', 'Mean_rec', 'Mean_uh'], 'Max':['Max', 'Max_rec','Max_uh'],
             'Min':['Min', 'Min_rec', 'Min_uh'] }

    assert len(labels)>=len(listdf)

    if norm=='h1':
        nm = f"$H^1$"
    else :
        nm = f"$L^2$"

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


def plotErrors(df, keys='Mean', norm='l2', texSave=False, name="plot.tex"):
    """plot nirb error given from a dataFrame keys take 'Min', 'Max' or 'Mean'.
        This will plot computed errors with and without rectification.
        The dataFrame df should be the result of function getDataStat()

    Args:
    -----
        df (pandas.dataFrame): _description_
        keys (str) : Min, Max or Mean. Defaults to Mean
        norm (str) : show which norm are displayed. Defaults to 'l2' norm
        texSave (bool) : say if save tex file or not. Defaults to False

    """

    xf = df.index

    key_list = {'Mean':['Mean', 'Mean_rec', 'Mean_uh'], 'Max':['Max', 'Max_rec','Max_uh'],
             'Min':['Min', 'Min_rec', 'Min_uh'] }

    normUHn = {'l2':r"$\Vert u^\mathcal{N}_h - u^N_{Hh}\Vert_{L^2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^N_{Hh}\Vert_{H^1}$"}
    normUhn = {'l2':r"$\Vert u^\mathcal{N}_h - u^N_{h}\Vert_{L^2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^N_{h}\Vert_{H^1}$"}
    normUh = {'l2':r"$\Vert u^\mathcal{N}_h - u^\mathcal{N}_{Hh}\Vert_{L^2}$" , 'h1': r"$\Vert u^\mathcal{N}_h - u^\mathcal{N}_{Hh}\Vert_{H^1}$"}

    keyUh = {'l2': 'l2(uh-uH)', 'h1':'h1(uh-uH)'}

    if norm=='h1':
        nm = f"$H^1$"
    else :
        nm = f"$L^2$"

    plt.scatter(xf, df[key_list[keys][0]], c='red', label=normUHn[norm] + 'w/o rect')
    plt.scatter(xf, df[key_list[keys][1]], c='blue', label=normUHn[norm] + 'w/ rect')
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
    if texSave :
        tikzplotlib.save(name)
        print(f"figure saved in {os.getcwd()}/{name}")
    plt.show()


def compare2dataFrame(df,dg, keys='Mean', norm='l2', dataLabel='pk', labellist=None, texSave=False):
    """Compare two dataFrames in respect of given keys and norm         keys take 'Min', 'Max' or 'Mean'
        the dataFrames df and dg should be the results of function getDataStat()

    Args:
    -----
        df (pandas.dataFrame): _description_
        dg (pandas.dataFrame): _description_
        keys (str) : Min, Max or Mean
        norm (str) : the norm to be plotted (l2 or h1)
        dataLabel (str) : the label of data to compare (pk : lagrange Pk interpolation, 'nested' : compare nested basis function
                        vs no nested one, 'greedy' : compare Greedy algo Vs no Greedy)
        texSave(bool) : Say if saving tex file or not. Defaults to False

    """

    xf = df.index
    xg = dg.index

    if norm=='h1' or norm=='H1':
        nm = f"$H^1$"
        normUHn = r"$\Vert u^\mathcal{N}_h - u^N_{Nh}\Vert_{H^1}$"
    else :
        nm = f"$L^2$"
        normUHn = r"$\Vert u^\mathcal{N}_h - u^N_{Nh}\Vert_{L^2}$"

    key_list = {'Mean':['Mean', 'Mean_rec', 'Mean_uh'], 'Max':['Max', 'Max_rec','Max_uh'],
             'Min':['Min', 'Min_rec', 'Min_uh'] }


    if labellist==None :
        if dataLabel=='greedy' :
            labeldf = ' w/ greedy'
            labeldg = ' w/o greedy'
        elif dataLabel=='nested':
            labeldf = ' Nested'
            labeldg = ' noNested'
        else :
            labeldf = '$\mathbb{P}_1$'
            labeldg = '$\mathbb{P}_2$'
    else :
        labeldf = labellist[0]
        labeldg = labellist[1]


    order = df[key_list[keys][0]]/dg[key_list[keys][0]]

    # plt.scatter(xf, df[key_list[keys][0]], marker='o', c='red', label=normUHn[norm] + ' w/o rect -'+ labeldf)
    # plt.scatter(xg, dg[key_list[keys][0]], marker='x', c='red', label=normUHn[norm] + ' w/o rect -'+ labeldg)

    # plt.scatter(xf, df[key_list[keys][1]], marker='o', c='blue', label=normUHn[norm] + ' w/ rectif -'+labeldf)
    # plt.scatter(xg, dg[key_list[keys][1]], marker='x', c='blue', label=normUHn[norm] + ' w/ rectif -'+labeldg)

    # plt.plot(xf, order, c='blue', label='Greedy/noGreedy')

    plt.plot(xf, df[key_list[keys][0]], c='red', label=labeldf + ' & w/o rect')
    plt.plot(xg, dg[key_list[keys][0]],'--', c='red', label=labeldg + ' & w/o rect')

    plt.plot(xf, df[key_list[keys][1]], c='blue', label=labeldf + ' & w/ rect' )
    plt.plot(xg, dg[key_list[keys][1]],'--', c='blue', label=labeldg + ' & w/ rect')

    plt.legend()
    plt.yscale('log')
    plt.xlabel("Number of basis function (N)")
    plt.ylabel(f"{nm} norm of Errors (in log scale)")
    plt.grid()
    if texSave :
        tikzplotlib.save(f"compareError{keys}.tex")
    # tikzplotlib.save(f"compareError{keys}.tex")
    plt.show()



if __name__ == "__main__":

    # path to csv file
    df_path = "errors20Params.csv"

    l2, h1 = getDataStat(pd.read_csv(df_path))

    l2.to_csv("l2.csv")
    h1.to_csv("h1.csv")

    # fig_l2 = l2.plot()
    # fig_l2.update_yaxes(type="log")
    # fig_l2.show()
    # tikzplotly.save("l2.tex", fig_l2)

    fig_h1 = h1.plot()
    # fig_h1.update_yaxes(type="log")
    # fig_h1.show()
    # tikzplotly.save("h1.tex", fig_h1)