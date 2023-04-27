import feelpp
import feelpp.mor as mor
import feelpp._core as core
import feelpp.toolboxes.core as tb_core
import feelpp.toolboxes.heat as heat
import sys, os
import pandas as pd
from feelpp.mor.online import Online
import subprocess
from numpy import abs



def loadParameterSpace(model_path):
    """load the parameter space from a model

    Parameters
    ----------
    model_path : str
        path to the model description (json file)

    Returns
    -------
    feelpp.mor._mor.ParameterSpace
        object containing the parameter space
    """
   
    crb_model_properties = mor.CRBModelProperties("", feelpp.Environment.worldCommPtr(), "")
    crb_model_properties.setup(model_path)
    Dmu = feelpp.mor._mor.ParameterSpace.New(crb_model_properties.parameters(), feelpp.Environment.worldCommPtr())
    return Dmu

def generate_sampling(model_path, Nsamples=1, samplingMode="random"):
    """Get sampling associated to the parameter space of the model

    Parameters
    ----------
    model_path : str
        path to the model description (json file)
    Nsamples : int, optional
        size of the sampling, by default 1
    samplingMode : str, optional
        samplingMode of sampling, by default "random"

    Returns
    -------
    feelpp.mor._mor.ParameterSpaceSampling
        sampling generated
    """
    Dmu = loadParameterSpace(model_path)
    s_heat = Dmu.sampling()
    s_heat.sampling(Nsamples, samplingMode)
    for i in range(Nsamples):
        # s_heat[i].setParameterNamed("h_bl", round(s_heat[i].parameterNamed("h_bl"), 2))
        # s_heat[i].setParameterNamed("h_amb", round(s_heat[i].parameterNamed("h_amb"), 2))
        s_heat[i].setParameterNamed("h_r", 6)
        # s_heat[i].setParameterNamed("T_bl", round(s_heat[i].parameterNamed("T_bl"), 2))
        # s_heat[i].setParameterNamed("T_amb", round(s_heat[i].parameterNamed("T_amb"), 2))
        # s_heat[i].setParameterNamed("E", round(s_heat[i].parameterNamed("E"), 2))
        # s_heat[i].setParameterNamed("k_lens", round(s_heat[i].parameterNamed("k_lens"), 2))
    return s_heat


def convert_model(mu, Dmu):
    mu_model = Dmu.element()
    k_lens = mu.parameterNamed("k_lens")
    h_amb = mu.parameterNamed("h_amb")
    h_bl = mu.parameterNamed("h_bl")
    h_r = mu.parameterNamed("h_r")
    T_amb = mu.parameterNamed("T_amb")
    T_bl = mu.parameterNamed("T_bl")
    E = mu.parameterNamed("E")

    mu_model.setParameterNamed("mu0", k_lens)
    mu_model.setParameterNamed("mu1", h_amb)
    mu_model.setParameterNamed("mu2", h_bl)
    mu_model.setParameterNamed("mu3", h_r)
  # mu_model.setParameterNamed("mu4", 1)
    mu_model.setParameterNamed("mu5", h_amb*T_amb + h_r*T_amb - E)
    mu_model.setParameterNamed("mu6", h_bl*T_bl)

    return mu_model


def convert_sampling(s_heat, Dmu):
    """Convert a sampling from the heat model to the eye2brain model

    Parameters
    ----------
    s_heat : feelpp.mor._mor.ParameterSpaceSampling
        Sampling of the heat model
    Dmu : feelpp.mor._mor.ParameterSpace
        Parameter space of the eye2brain model

    Returns
    -------
    feelpp.mor._mor.ParameterSpaceSampling
        Sampling of the eye2brain model
    """
    s = Dmu.sampling()
    N = len(s_heat)
    s.sampling(N, "random")
    for i in range(N):
        mu = s_heat[i]
        mu_model = convert_model(mu, Dmu)
        s[i] = mu_model
    return s


def convert_results_to_dataframe(res):
    """Convert the result of the online phase to dataframe

    Parameters
    ----------
    res : list of feelpp.mor._mor.CRBResults
        list of results

    Returns
    -------
    pandas.DataFrame
        dataframe with the results : parameter, output, errorBound
    """
    names = res[0].parameter().parameterNames()

    df = pd.DataFrame(columns = names + ["RB_output", "errorBound"])
    for i, r in enumerate(res):
        mu = r.parameter()
        df.loc[i] = [mu.parameterNamed(n) for n in names] + [r.output(), r.errorBound()]
    return df

def convert_parameters_to_dataframe(inputs):
    """Convert the result of the online phase to dataframe

    Parameters
    ----------
    inputs : list of feelpp.mor._mor.ParameterSpaceSampling
        list of inputs from the heat model

    Returns
    -------
    pandas.DataFrame
        dataframe with the results : parameter
    """
    names_heat = inputs[0].parameterNames()

    df = pd.DataFrame(columns = names_heat)
    for i in range(len(inputs.getVector())):
        mu_heat = inputs[i]
        df.loc[i] = [mu_heat.parameterNamed(n) for n in names_heat]
    return df

############################################################################################################
# PFEM part

def run_pfem(sample):
    if feelpp.Environment.isMasterRank():
        print("Offline phase")
        np = feelpp.Environment.numberOfProcessors()
        APP_DIR = "~/Documents/code/feelpp-dev/build/mor/mor/examples/eye2brain/"
        app_path = os.path.join(APP_DIR, "feelpp_mor_eye2brainapp")
        cfg_path = os.path.join(os.path.join( os.path.dirname(os.path.abspath(__file__)), "eye2brain/eye2brain.cfg"))
        out = []
        for i, mu in enumerate(sample.getVector()):
            mu_str = " ".join([f"{mu.parameterNamed(n)} " for n in mu.parameterNames()])[:-1]
            command = f"{app_path} --config-file {cfg_path} --eye2brain.run.mode 0 --crb.user-parameters \"{mu_str}\""
            print(f"Running command : {command}")
            output = subprocess.getoutput(command)
            log_path = f"logs/output_{i}.log"
            f = open(log_path, "w")
            f.write(output)
            f.close()
            o = subprocess.getoutput(f"cat {log_path} | grep Eye2Brain::output")
            out.append( float(o.split(' ')[-1]) )
    return out


############################################################################################################
# FEM simulation using the toolbox

def assembleToolbox(tb, mu):
    """Assemble the toolbox with the given parameters

    Parameters
    ----------
    tb : feelpp.toolboxes._toolboxes.Toolbox
        toolbox to assemble
    mu : feelpp.mor._mor.Parameter
        parameters to use
    """
    for i in range(0,mu.size()):
        tb.addParameterInModelProperties(mu.parameterName(i), mu(i))

    for i in range(0,mu.size()):
        tb.addParameterInModelProperties(mu.parameterName(i), mu(i))

    tb.updateParameterValues()

def run_toolbox(app, sample, output_key = "Statistics_cornea_mean"):
    CFG_DIR = os.path.join( os.path.dirname(os.path.abspath(__file__)), "eye2brain" )
    cfg_file = os.path.join(CFG_DIR, "eye-linear.cfg")

    app.setConfigFile(cfg_file)

    heatBox = heat.heat(dim=3, order=1)
    heatBox.init()

    out = []
    for i, mu in enumerate(sample.getVector()):
        assembleToolbox(heatBox, mu)
        heatBox.solve()
        heatBox.exportResults()

        l = heatBox.postProcessMeasures()
        meas = heatBox.postProcessMeasures().values()

        if feelpp.Environment.isMasterRank():
            print(i, "mu meas =", meas[output_key])
        out.append(meas[output_key])
    
    return out



if __name__ == '__main__':

    opts = tb_core.toolboxes_options("heat").add(mor.makeToolboxMorOptions()).add(mor.makeCRBOptions())
    config = feelpp.globalRepository(f"eye2brain/compFE-RB")
    app = feelpp.Environment(sys.argv, opts=opts, config=config)
    crbdir = app.repository().globalRoot().string()
    m_def = os.path.join( os.path.dirname(os.path.abspath(__file__)), "eye2brain/crb_param.json" )
    name = "eye2brain_p1g1"

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-N', '--Num', dest='N', help='number of parameters to evaluate the plugin', type=int, default=1)
    parser.add_option('-d', '--dir', dest='dir', help='directory location of crbdb directory', default=crbdir)
    parser.add_option('-n', '--name', dest='name', help='name of the plugin', type="string", default=name)
    parser.add_option('-i', '--id', dest='dbid', help='DB id to be used', type="string")
    parser.add_option('-l', '--load', dest='load', help='type of data to be loadedoaded', type="string",default="rb")
    parser.add_option('-m', '--model', dest='model', help='path to the model description', type="string", default=m_def)

    parser.add_option('-t', '--toolbox-only', dest='toolbox', help='run only toolbox output for given sampling', default=False)
    parser.add_option('-s', '--sampling', dest='sampling', help='path to the heat model sampling file (if None, one is generated)', type="string", default=None)
    parser.add_option('-o', '--output', dest='output', help='path to the output file', type="string", default=None)

    (options, args) = parser.parse_args()


    if options.sampling is None:
        s_heat = generate_sampling(options.model, 0)
    else:
        Dmu = loadParameterSpace(options.model)
        s_heat = Dmu.sampling()
        N_read = s_heat.readFromFile(options.sampling)
        if feelpp.Environment.isMasterRank():
            print("Read", N_read, "samples from", options.sampling)
    
    df = convert_parameters_to_dataframe(s_heat)

    if not options.toolbox:
        o = Online(options.name, options.dir)
        s = convert_sampling(s_heat, o.rbmodel.parameterSpace())

        res_rb = o.run(s)
        df_rb = convert_results_to_dataframe(res_rb)

        df = pd.concat([df, df_rb], axis=1)

        if not feelpp.Environment.isParallel():
            out = run_pfem(s)
            df['PFEM_output'] = out
        else:
            if feelpp.Environment.isMasterRank():
                print("Cannot run parallel of PFEM, sorry !")

    else:
        if options.output is not None:
            df0 = pd.read_csv(options.output)
        else:
            df0 = None

        df = pd.concat([df, df0], axis=1)


    out = run_toolbox(app, s_heat)
    df['FEM_output'] = out

    if feelpp.Environment.isMasterRank():
        print(df)
        print()
        if 'errorBound' in df.columns:
            print("Error bounds")
            print(df['errorBound'])
        if 'RB_output' in df.columns:
            print("\nRB - PFEM")
            print(abs(df["RB_output"] - df["FEM_output"]))
        if 'PFEM_output' in df.columns:
            print("\nFEM - PFEM")
            print(abs(df["FEM_output"] - df["PFEM_output"]))

    df.to_csv("results.csv")

    sys.exit(0)