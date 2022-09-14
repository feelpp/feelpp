from nirb import *

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

    doRectification=True
    computeCoarse=False 
    if doRectification:
        computeCoarse=True 

    nirb_model = NIRB(dim, H, h, toolboxesOptions, cfg_path, model_path, geo_path, doRectification=doRectification)

    nirb_model.initProblem(3,computeCoarse=computeCoarse)
    nirb_model.generateOperators()
    nirb_model.generateReducedBasis()

    # nirb_model.orthonormalizeL2()
    # nirb_model.orthonormalizeH1()
    nirb_model.saveData()
    print("Is L2 orthonormalized ?", nirb_model.checkL2Orthonormalized())
    print("Is H1 orthonormalized ? ", nirb_model.checkH1Orthonormalized())

    print(f"[NIRB] Offline part Done !")

    nirb_model.loadData()
    uc = nirb_model.getInterpSol()
    uc = nirb_model.getCompressedSol()
    nirb_model.getOnlineSol() 
    print("[nirb main] online sol befor error", nirb_model.onlineSol.to_petsc().vec()[:])
    error = nirb_model.computeErrors() 
    print("[nirb main] online sol after error", nirb_model.onlineSol.to_petsc().vec()[:])

    print("----------- Errors ----------------------")
    print(f"Nb mode = {error[0]}")
    print(f"L2 norm= {error[1]}")
    print(f"Inf norm = {error[2]}")

    print(f"[NIRB] Online part Done !!")