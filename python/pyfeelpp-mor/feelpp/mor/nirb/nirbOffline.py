from asyncio.base_futures import _FINISHED
from tracemalloc import start
from turtle import st
from nirb import *
from utils import WriteVecAppend 
import timeit 
from timeit import default_timer as timer
import time 

if __name__ == "__main__":

    # fineness of two grids
    H = 0.1  # CoarseMeshSize
    h = H**2 # Fine mesh size
    dim = 2

    PWD = os.getcwd()
    toolboxesOptions='heat'
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}{modelfile[toolboxesOptions]}.cfg"
    geo_path = f"{modelsFolder}{modelfile[toolboxesOptions]}.geo"
    model_path = f"{modelsFolder}{modelfile[toolboxesOptions]}.json"

    doRectification=True 
    nbSnap = 10
    if len(sys.argv)>1:
        nbSnap = int(sys.argv[1])

    
    # star=timeit.timeit()
    star = time.perf_counter()

    nirb_off = nirbOffline(dim, H, h, toolboxesOptions, cfg_path, model_path, geo_path, doRectification=doRectification)

    nirb_off.initProblem(nbSnap)
    nirb_off.generateOperators()
    nirb_off.generateReducedBasis()

    # nirb_off.BiOrthonormalization()
    
    # nirb_model.orthonormalizeL2()
    # nirb_model.orthonormalizeH1()
    nirb_off.saveData()
    print("Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    print("Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())


    finish = time.perf_counter()

    perf = []
    perf.append(nbSnap)
    perf.append(finish-star)

    file='nirbOffline_time_exec.txt'
    WriteVecAppend(file,perf)

    print(f"[NIRB] Offline Elapsed time = ", finish-star)
    print(f"[NIRB] Offline part Done !")
