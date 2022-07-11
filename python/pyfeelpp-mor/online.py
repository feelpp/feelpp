import sys, os
import feelpp.mor.reducedbasis.reducedbasis as mor_rb
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp

HOME = os.environ['HOME']

if __name__ == '__main__':
    config_file = f'{HOME}/feel/feelpp/python/pyfeelpp-mor/tests/testcase/square/2d/testcase2d.cfg'
    odir = 'testcase2d'
    dim = 2
    time_dependant = False

    split = config_file.split('/')
    case = '/'.join(split[:-1])
    dir = odir
    casefile = split[-1]

    config = feelpp.globalRepository("generate_basis-testcase2d")
    sys.argv = ['online']
    o = toolboxes_options("heat")
    o.add(makeToolboxMorOptions())

    e = feelpp.Environment(sys.argv, opts=o, config=config)

    feelpp.Environment.setConfigFile(f'{case}/{casefile}')

    print("current working directory was", os.getcwd())
    os.chdir(f"{HOME}/feel/feelpp/generate_basis-testcase2d/np_1/")
    print("current working directory is now", os.getcwd())

    name =  "generate_basis-testcase2d"# + "-np_" +  str(feelpp.Environment.numberOfProcessors())
    name = name.replace('--', '-')

    # Set the toolboxes
    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    # model = toolboxmor_2d() if DIM==2 else toolboxmor_3d()
    model = toolboxmor(name=name, dim=dim, time_dependent=time_dependant)
    model.setFunctionSpaces( Vh=heatBox.spaceTemperature() )


    def assembleDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleRhs()

    def assembleMDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleMatrix()

    model.setAssembleDEIM(fct=assembleDEIM)

    model.setAssembleDEIM(fct=assembleDEIM)
    model.setAssembleMDEIM(fct=assembleMDEIM)

    model.initModel()

    heatBoxDEIM = heat(dim=dim, order=1)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    heatBoxMDEIM = heat(dim=dim, order=1)
    meshMDEIM = model.getMDEIMReducedMesh()
    heatBoxMDEIM.setMesh(meshMDEIM)
    heatBoxMDEIM.init()

    def assembleOnlineDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxDEIM.updateParameterValues()
        return heatBoxDEIM.assembleRhs()

    def assembleOnlineMDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxMDEIM.updateParameterValues()
        return heatBoxMDEIM.assembleMatrix()

    model.setOnlineAssembleDEIM(assembleOnlineDEIM)
    model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

    model.postInitModel()
    model.setInitialized(True)

    # Run online computation
    Dmu = model.parameterSpace()

    names = Dmu.parameterNames()
    mumax = Dmu.mumax()
    mumin = Dmu.mumin()

    print("       ", end='')
    for name in names:
        print(f"{name:<9}", end="")
    print()
    print("mumin", mumin)
    print("mumax", mumax)

    # Here, we create a dict containing the bounds for each parameter
    def ints_of_model(Dmu):
        mumax = Dmu.mumax()
        mumin = Dmu.mumin()
        names = Dmu.parameterNames()
        ints = {}
        for name in names:
            ints[name] = {'min': mumin.parameterNamed(name), 'max': mumax.parameterNamed(name)}
        return ints
    ints = ints_of_model(Dmu)
    ints

    basis = mor_rb.reducedbasis(None)

    basis.loadReducedBasis(f'{HOME}/feel/feelpp/generate_basis-testcase2d/np_1/reducedbasis.json', model)

    mu = basis.mubar
    mu.view()

    print("basis.N_output : ", basis.N_output)
    print("basis.output_names : ", basis.output_names)

    def f(k, k_2):
        mu = Dmu.element(True,False)
        mu.setParameters({"k_2":k_2})
        print("[k_2]")
        print(mu)

        uN,sN = basis.getSolutions(mu, k=k)
        print(f"Output : ",basis.output_names[k], sN)

    f(0,10)
    f(1,10)
    f(2,10)
    f(3,10)
