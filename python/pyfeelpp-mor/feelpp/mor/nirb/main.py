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

    nirb_model = NIRB(dim, H, h, toolboxesOptions, cfg_path, model_path, geo_path, doRectification=False)

    nirb_model.initProblem(10)
    nirb_model.generateOperators()
    nirb_model.generateReducedBasis()

    # nirb_model.saveData()
    nirb_model.checkOrthonormalized()