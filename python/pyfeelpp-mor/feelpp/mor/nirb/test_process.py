#%%
## Import module 
from postProcess import * 

#%%

if __name__ == "__main__":
    
    import sys 

    dataPath = "/data/home/elarif/feelppdb/nirb/heat/np_1/results/Rect/"

    #%% 
    # Get dataFrame from csv file 
    # file = dataPath + "errors50ParamsLambda10P1.csv"
    file = dataPath + "noGreedy/errors50Params.csv"
    file = dataPath + "noGreedy/errors30Params_s4.csv"
    # file = dataPath + "noGreedy/offlineError.csv"
    
    # load absolute errors
    # file = dataPath + "noGreedy/offlineConvError.csv"
    dfGlob = pd.read_csv(file, sep=',')
    
    # # compute relative errors 
    # dfRel = getRelativeErrors(dfGlob)

    ### Get stats for all parameters 
    l2df  = getDataStat(dfGlob, h1norm=False) # l1 and h1 error associated 
    # print(l2df.head())
    # l2dfRel, h1dfRel = getDataStat(dfRel) # // 

    plotErrors(l2df, keys='Mean')

    # %%
    ## Tronctae error datas into N parameter (N<=50)
    # N = 5
    # dfN = troncateNparam(dfGlob, N) 
    
    # l2dfN, h1dfN   = getDataStat(dfN) # l1 and h1 error associated with statistical infos 
    
    # plot_dataFrame(l2dfN)

    #%%
    ### Compare list of error files according to (\lambda, parameters, ...) 
    listlabel = []
    listdfl2, listdfh1 = [], []
    indexes = [2, 4, 9, 16, 25]
    indexes = [0, 1, 3, 5, 7, 9, 10, 11, 12, 17]
    sk = 's9' # s9
    for st in indexes:
        file = dataPath + f"noGreedy/errors30Params_s4lmd{st}.csv"
        dfG = pd.read_csv(file, sep=',')
        dl2, dh1 = getDataStat(dfG)
        listdfl2.append(dl2)
        listdfh1.append(dh1)
        listlabel.append(f"$\lambda = 1.e({-st})$")
        # listlabel.append(f"{st} param")
    
    compareListOfDataFrams(listdfl2, norm='Min', listlabel=listlabel, rectif=True)

    # %%
    # Compare two data frames (P1 vs P2) or (Greedy vs noGreedy)  
    
    # dataPath = "/data/home/elarif/feelppdb/nirb/heat/np_1/results/Rect/"

    filep1 = dataPath + "noGreedy/errors10ParamsEmboited.csv"
    filep2 = dataPath + "noGreedy/errors10ParamsNoEmboited.csv"

    filep1 = dataPath + "Greedy/errors20Params.csv"
    filep2 = dataPath + "noGreedy/errors20Params.csv"

    filep1 = dataPath + "noGreedy/errors10ParamsLU.csv"
    filep2 = dataPath + "noGreedy/errors10ParamsGamg.csv"

    dfp1 = pd.read_csv(filep1, sep=',')
    dfRelp1 = getRelativeErrors(dfp1)

    dfp2 = pd.read_csv(filep2, sep=',')
    dfRelp2 = getRelativeErrors(dfp2)

    l2p1, h1p1 = getDataStat(dfp1)
    l2p2, h1p2 = getDataStat(dfp2)

    labellist = ['lu', 'gamg']
    compare2dataFrame(l2p1, l2p2, dataLabel='greedy', keys='Max', labellist=labellist)

#%%