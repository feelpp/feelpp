# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib

# %%

def plot_error(dirs, names, Ns):
    """Plot errors obtained from many directories

    Args:
        dirs (list of str): list of path to export directories
        names (list of str): list of names for the legend
        Ns (list of int): sizes of the basis used for the computation of the errors
    """
    D = len(dirs)
    assert(D == len(names))
    NN = len(Ns)

    cs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

    fig, ax = plt.subplots(1, 2, figsize=(20,10))

    errIntL2 = 0
    errIntLinf = 0

    for d in range(D):
        minsL2  = np.zeros(NN)
        meansL2 = np.zeros(NN)
        maxsL2  = np.zeros(NN)

        minsLinf  = np.zeros(NN)
        meansLinf = np.zeros(NN)
        maxsLinf  = np.zeros(NN)

        for i, N in enumerate(Ns):
            df = pd.read_csv(f"{dirs[d]}/errorParams/errors{N}.csv")
            minsL2[i]  = df["l2u-uHn"].min()
            meansL2[i] = df["l2u-uHn"].mean()
            maxsL2[i]  = df["l2u-uHn"].max()

            minsLinf[i]  = df["lINFu-uHn"].min()
            meansLinf[i] = df["lINFu-uHn"].mean()
            maxsLinf[i]  = df["lINFu-uHn"].max()

            errIntL2 += df["l2u-uH"].mean()
            errIntLinf += df["lINFu-uH"].mean()

            print(f"[NIRB online] {names[d]}: N={N},\n\tmin L2 error = {minsL2[i]},\n\tmean L2 error = {meansL2[i]},\n\tmax L2 error = {maxsL2[i]}")
            print(f"[NIRB online] {names[d]}: N={N},\n\tmin Linf error = {minsLinf[i]},\n\tmean Linf error = {meansLinf[i]},\n\tmax Linf error = {maxsLinf[i]}")

        ax[0].plot(Ns, minsL2, c=cs[d], label=f"L2 {names[d]}")
        ax[0].plot(Ns, meansL2, '--', c=cs[d])
        ax[0].plot(Ns, maxsL2, c=cs[d])
        ax[0].fill_between(Ns, minsL2, maxsL2, alpha=0.2, color=cs[d])

        ax[1].plot(Ns, minsLinf, c=cs[d], label=f"Linf {names[d]}")
        ax[1].plot(Ns, meansLinf, '--', c=cs[d])
        ax[1].plot(Ns, maxsLinf, c=cs[d])
        ax[1].fill_between(Ns, minsLinf, maxsLinf, alpha=0.2, color=cs[d])

    ax[0].axhline(errIntL2 / D, c='k', linestyle='--', label=r"$\Vert u_H - u_h\Vert$")
    ax[1].axhline(errIntLinf / D, c='k', linestyle='--', label=r"$\Vert u_H - u_h\Vert$")

    # Plot
    ax[0].set_title(r"Error $L^2$")
    ax[0].set_xlabel("N")
    ax[0].set_ylabel("Error")
    # ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(Ns[0], Ns[-1])
    ax[0].legend()

    ax[1].set_title(r"Error $L^\infty$")
    ax[1].set_xlabel("N")
    ax[1].set_ylabel("Error")
    # ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(Ns[0], Ns[-1])
    ax[1].legend()

    # Save
    # tikzplotlib.save("plot.tex")
    plt.show()



# %%
def plot_time(csv_file):
    df = pd.read_csv(csv_file)
    plt.axhline(df['time_toolbox'].mean(), c='k', linestyle='--', label='Fine toolbox')
    plt.plot(df['N'], df['time_nirb'], '+-', label='NIRB w/o rect.')
    plt.plot(df['N'], df['time_nirb_rect'], '+-', label='NIRB w/ rect.')

    plt.xlabel('N')
    plt.ylabel('Time (s)')
    # plt.xscale('log')
    plt.yscale('log')

    plt.legend()
    plt.show()
    # tikzplotlib.save("plot.tex")

# %%
if __name__ == "__main__":
    import sys
    plot_error(sys.argv[1], sys.argv[2:])
    plot_time(sys.argv[1])
