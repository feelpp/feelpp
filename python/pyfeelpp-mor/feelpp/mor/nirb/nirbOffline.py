import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import time 

if __name__ == "__main__":

    # fineness of two grids
    H = 0.25  # CoarseMeshSize
    h = H**2 # Fine mesh size
    dim = 2

    PWD = os.getcwd()
    toolboxType='heat'
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}{modelfile[toolboxType]}.cfg"
    geo_path = f"{modelsFolder}{modelfile[toolboxType]}.geo"
    model_path = f"{modelsFolder}{modelfile[toolboxType]}.json"

    e = init_feelpp_environment(toolboxType, cfg_path)

    doRectification=True  
    nbSnap = 10
    if len(sys.argv)>=2:
        nbSnap = int(sys.argv[1])
        if len(sys.argv)>=3: 
            H = float(sys.argv[2])
            h = H**2

    star=time.time()

    nirb_off = nirbOffline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, doRectification=doRectification)

    nirb_off.initProblem(nbSnap)
    nirb_off.generateOperators()
    nirb_off.generateReducedBasis(regulParam=1.e-10)
    # nirb_off.BiOrthonormalization()

    # nirb_off.orthonormalizeL2()
    # nirb_off.orthonormalizeH1()
    nirb_off.saveData()

    
    print("Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    print("Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())

    # finish = time()
    finish = time.time()
    


    perf = []
    perf.append(nbSnap)
    perf.append(finish-star)

    file='nirbOffline_time_exec.txt'
    WriteVecAppend(file,perf)

    print(f"[NIRB] Offline Elapsed time = ", finish-star)
    print(f"[NIRB] Offline part Done !")
