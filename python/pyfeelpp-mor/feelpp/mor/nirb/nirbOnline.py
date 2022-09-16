from timeit import timeit
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend

if __name__ == "__main__":

    # fineness of two grids
    H = 0.1  # CoarseMeshSize
    h = H**2 # fine mesh size
    dim = 2

    PWD = os.getcwd()
    toolboxesOptions='heat'
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}{modelfile[toolboxesOptions]}.cfg"
    geo_path = f"{modelsFolder}{modelfile[toolboxesOptions]}.geo"
    model_path = f"{modelsFolder}{modelfile[toolboxesOptions]}.json"

    start=timeit() 
    doRectification=False

    nirb_on = nirbOnline(dim, H, h, toolboxesOptions, cfg_path, model_path, geo_path, doRectification=doRectification)

    nirb_on.loadData()
    nirb_on.getInterpSol()
    nirb_on.getCompressedSol()
    nirb_on.getOnlineSol()
    online1 = nirb_on.onlineSol.to_petsc().vec()[:] # en commentant cette ligne ça produite des nan à la solution onlineSol après computeErrors

    error = nirb_on.computeErrors()

    online2 = nirb_on.onlineSol.to_petsc().vec()[:]

    print("[nirb main] diff online sol befor/after error", np.max(np.abs(online1-online2)))

    print("----------- Errors ----------------------")
    print(f"Nb mode = {error[0]}")
    print(f"L2 norm= {error[1]}")
    print(f"Inf norm = {error[2]}")


    finish = timeit() 

    perf = []
    perf.append(nirb_on.N)
    perf.append(finish-start)

    file='nirbOnline_time_exec.txt'
    WriteVecAppend(file,perf)

    if doRectification:
        file='nirb_error_rectif.txt'
    else :
        file='nirb_error.txt'
    WriteVecAppend(file,error)

    print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")