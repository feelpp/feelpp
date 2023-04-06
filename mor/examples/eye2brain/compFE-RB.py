import feelpp
import feelpp.mor as mor
import sys, os
import pandas as pd
from feelpp.mor.online import Online
import subprocess


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
    return s_heat


def convert_model(mu, Dmu):
    mu_model = Dmu.element()
    k_lens = mu.parameterNamed("k_lens")
    h_amb = mu.parameterNamed("h_amb")
    h_bl = mu.parameterNamed("h_bl")
    T_amb = mu.parameterNamed("T_amb")
    T_bl = mu.parameterNamed("T_bl")
    E = mu.parameterNamed("E")

    mu_model.setParameterNamed("mu0", k_lens)
    mu_model.setParameterNamed("mu1", h_amb)
    mu_model.setParameterNamed("mu2", h_bl)
  # mu_model.setParameterNamed("mu3", 1)
    mu_model.setParameterNamed("mu4", h_amb*T_amb + 6*T_amb + E)
    mu_model.setParameterNamed("mu5", h_bl*T_bl)
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


def convert_to_dataframe(res, inputs):
    """Convert the result of the online phase to dataframe

    Parameters
    ----------
    res : list of feelpp.mor._mor.CRBResults
        list of results
    inputs : list of feelpp.mor._mor.ParameterSpaceSampling
        list of inputs from the heat model

    Returns
    -------
    pandas.DataFrame
        dataframe with the results : parameter, output, errorBound
    """
    names_heat = inputs[0].parameterNames()
    names = res[0].parameter().parameterNames()

    df = pd.DataFrame(columns=names_heat + names + ["RB_output", "errorBound"])
    for i, r in enumerate(res):
        mu_heat = inputs[i]
        mu = r.parameter()
        df.loc[i] = [mu_heat.parameterNamed(n) for n in names_heat] + [mu.parameterNamed(n) for n in names] + [r.output(), r.errorBound()]
    return df


def run_offline(sample):
    APP_DIR = "~/Documents/code/feelpp-dev/build/mor/mor/examples/eye2brain/"
    app_path = os.path.join(APP_DIR, "feelpp_mor_eye2brainapp")
    cfg_path = os.path.join(os.path.join( os.path.dirname(os.path.abspath(__file__)), "eye2brain/eye2brain.cfg"))
    out = []
    for i, mu in enumerate(sample.getVector()):
        mu_str = " ".join([f"{mu.parameterNamed(n)} " for n in mu.parameterNames()])[:-1]
        command = f"{app_path} --config-file {cfg_path} --eye2brain.run.mode 0 --crb.user-parameters \"{mu_str}\""
        output = subprocess.getoutput(command)
        log_path = f"output_{i}.log"
        f = open(log_path, "w")
        f.write(output)
        f.close()
        o = subprocess.getoutput(f"cat {log_path} | grep Eye2Brain::output")
        out.append( float(o.split(' ')[-1]) )
    return out


if __name__ == '__main__':

    app = feelpp.Environment(sys.argv, opts=mor.makeCRBOptions())
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
    
    (options, args) = parser.parse_args()
    print("members:", mor.CRBLoad.__members__)
    print("rb members:", mor.CRBLoad.__members__["rb"])
    print("crbdir:", crbdir)

    o = Online(options.name, options.dir)

    s_heat = generate_sampling(options.model, options.N)
    s = convert_sampling(s_heat, o.rbmodel.parameterSpace())

    res_rb = o.run(s)

    df = convert_to_dataframe(res_rb, s_heat)
    print(df)

    out = run_offline(s)
    df['PFEM_output'] = out

    print(df)
    df.to_csv("results.csv")

    sys.exit(0)