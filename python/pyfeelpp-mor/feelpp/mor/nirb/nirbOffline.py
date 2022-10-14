import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import time
import json

if __name__ == "__main__":

    dim = 2
    if dim == 2:
        H = 0.1  # CoarseMeshSize
        h = H**2 # Fine mesh size
    else:
        # fineness of two grids
        H = 0.5  # CoarseMeshSize
        h = 0.1  # Fine mesh size

    PWD = os.getcwd()
    toolboxType='heat'
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}{modelfile[toolboxType]}.cfg"
    geo_path = f"{modelsFolder}{modelfile[toolboxType]}.geo"
    model_path = f"{modelsFolder}{modelfile[toolboxType]}.json"
    if dim == 3:
        cfg_path = f"{modelsFolder}/thermal-fin-3d/thermal-fin.cfg"
        geo_path = f"{modelsFolder}/thermal-fin-3d/fin.geo"
        model_path = f"{modelsFolder}/thermal-fin-3d/thermal-fin.json"

    e = init_feelpp_environment(toolboxType, cfg_path)

    doRectification=True
    nbSnap = 5
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
    nirb_off.saveData(force=True)


    # ######################################################################
    # ### test clone 
    # b = feelpp.backend(worldcomm=feelpp.Environment.worldCommPtr())
    # v = b.newVector(dm=nirb_off.Xh.mapPtr())

    # v.setConstant(1.0)


    # print("------------------------------------------")
    # print('type v :', type(v))
    # vv = v.clone()
    # print('type v clone ', type(vv))
    # print('val v', v.vec()[:])
    # print('val v clone', vv.vec()[:])


    # print("------------------------------------------")
    # w = nirb_off.Xh.element(v)
    # print('type w:', type(w))
    # w.setConstant(15.)
    # ww = w.clone()
    # print('type w clone ', type(ww))
    # print('val w', w.to_petsc().vec()[:])
    # print('val w clone', ww.to_petsc().vec()[:])

    # print("------------------------------------------")

    # x = w.clone()
    # print('type x:', type(x))
    # x.setConstant(30.)
    # xx = x.clone()
    # print('type x clone ', type(xx))
    # print('val x', x.to_petsc().vec()[:])
    # print('val x clone', xx.to_petsc().vec()[:])
    # print("------------------------------------------")

    # ######################################################################
    # ### test addVector 
    # b = feelpp.backend(worldcomm=feelpp.Environment.worldCommPtr())
    # vv = b.newVector(dm=nirb_off.Xh.mapPtr()) # a petsc vector 
    # ww = nirb_off.Xh.element(vv) # a Pch element
    # kk = ww.clone()              # an Ublas vector 
    
    # v = vv

    # u = nirb_off.reducedBasis[0]

    # print('type u', type(u))
    # print('type v :', type(v))

    # mat = nirb_off.l2ScalarProductMatrix
    # mat.zero()
    # ww.setConstant(1.0)
    # mat.setDiagonal(ww)
    # print('mat =', mat.mat()[:,:])

    # v.setConstant(1.0)
    # u.setConstant(2.0)

    # print("v before addVector :", v.to_petsc().vec()[:])
    # print('u before addVector :', u.to_petsc().vec()[:])

    # print("------------------------------------------")
    # v.addVector(u, mat)
    # print('u after addVector', u.to_petsc().vec()[:])
    # print('v after addVector', v.to_petsc().vec()[:])
    # print('type v after addVector', type(v))
    # print("------------------------------------------")
    
    # ######################################################################
    
    finish = time.time()

    # nirb_off.orthonormalizeL2()
    # nirb_off.orthonormalizeH1()
    # nirb_off.saveData()

    # print("Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    # print("Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())

    perf = []
    perf.append(nbSnap)
    perf.append(finish-star)

    info = nirb_off.getInformations()

    if feelpp.Environment.isMasterRank():
        print(json.dumps(info, sort_keys=True, indent=4))
        print(f"[NIRB] Offline Elapsed time = ", finish-star)
        print(f"[NIRB] Offline part Done !")
