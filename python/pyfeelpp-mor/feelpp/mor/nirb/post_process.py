# %%
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib

# %%

def plot_error(dir_name, Ns):
    N = len(Ns)
    # l2u-uHn, lINFu-uHn, l2u-uH, lINFu-uH
    mins2Nirb_norect = np.zeros(N)
    maxs2Nirb_norect = np.zeros(N)
    moy2Nirb_norect  = np.zeros(N)

    minsInfNirb_norect = np.zeros(N)
    maxsInfNirb_norect = np.zeros(N)
    moyInfNirb_norect  = np.zeros(N)

    mins2Int_norect = np.zeros(N)
    maxs2Int_norect = np.zeros(N)
    moy2Int_norect  = np.zeros(N)

    minsInfInt_norect = np.zeros(N)
    maxsInfInt_norect = np.zeros(N)
    moyInfInt_norect  = np.zeros(N)


    mins2Nirb_rect = np.zeros(N)
    maxs2Nirb_rect = np.zeros(N)
    moy2Nirb_rect  = np.zeros(N)

    minsInfNirb_rect = np.zeros(N)
    maxsInfNirb_rect = np.zeros(N)
    moyInfNirb_rect  = np.zeros(N)

    mins2Int_rect = np.zeros(N)
    maxs2Int_rect = np.zeros(N)
    moy2Int_rect  = np.zeros(N)

    minsInfInt_rect = np.zeros(N)
    maxsInfInt_rect = np.zeros(N)
    moyInfInt_rect  = np.zeros(N)

    for i,N in enumerate(Ns):
        df_norect = pd.read_csv(f"{dir_name}/errorParams/errors{N}.csv")
        mins2Nirb_norect[i] = df_norect["l2u-uHn"].min()
        maxs2Nirb_norect[i] = df_norect["l2u-uHn"].max()
        moy2Nirb_norect[i]  = df_norect["l2u-uHn"].mean()

        minsInfNirb_norect[i] = df_norect["lINFu-uHn"].min()
        maxsInfNirb_norect[i] = df_norect["lINFu-uHn"].max()
        moyInfNirb_norect[i]  = df_norect["lINFu-uHn"].mean()

        mins2Int_norect[i] = df_norect["l2u-uH"].min()
        maxs2Int_norect[i] = df_norect["l2u-uH"].max()
        moy2Int_norect[i]  = df_norect["l2u-uH"].mean()

        minsInfInt_norect[i] = df_norect["lINFu-uH"].min()
        maxsInfInt_norect[i] = df_norect["lINFu-uH"].max()
        moyInfInt_norect[i]  = df_norect["lINFu-uH"].mean()

        df_rect = pd.read_csv(f"{dir_name}/errorParamsRectif/errors{N}.csv")
        mins2Nirb_rect[i] = df_rect["l2u-uHn"].min()
        maxs2Nirb_rect[i] = df_rect["l2u-uHn"].max()
        moy2Nirb_rect[i]  = df_rect["l2u-uHn"].mean()

        minsInfNirb_rect[i] = df_rect["lINFu-uHn"].min()
        maxsInfNirb_rect[i] = df_rect["lINFu-uHn"].max()
        moyInfNirb_rect[i]  = df_rect["lINFu-uHn"].mean()

        mins2Int_rect[i] = df_rect["l2u-uH"].min()
        maxs2Int_rect[i] = df_rect["l2u-uH"].max()
        moy2Int_rect[i]  = df_rect["l2u-uH"].mean()

        minsInfInt_rect[i] = df_rect["lINFu-uH"].min()
        maxsInfInt_rect[i] = df_rect["lINFu-uH"].max()
        moyInfInt_rect[i]  = df_rect["lINFu-uH"].mean()

    print("MinNitb_norect", mins2Nirb_norect)
    print("MaxNitb_norect", maxs2Nirb_norect)
    print("MoyNitb_norect", moy2Nirb_norect)

    print("MinInfNitb_norect", minsInfNirb_norect)
    print("MaxInfNitb_norect", maxsInfNirb_norect)
    print("MoyInfNitb_norect", moyInfNirb_norect)

    print("MinInt_rect", mins2Int_rect)
    print("MaxInt_rect", maxs2Int_rect)
    print("MoyInt_rect", moy2Int_rect)

    print("MinInfInt_rect", minsInfInt_rect)
    print("MaxInfInt_rect", maxsInfInt_rect)
    print("MoyInfInt_rect", moyInfInt_rect)


    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(20,10))

    ax[0].plot(Ns, mins2Nirb_norect, c="tab:blue", label=r"$\Vert u_h - u_{hH}^N\Vert_{L^2}$ w/o rect.")
    ax[0].plot(Ns, maxs2Nirb_norect, c="tab:blue")
    ax[0].plot(Ns, moy2Nirb_norect, c="tab:blue")
    ax[0].fill_between(Ns, mins2Nirb_norect, maxs2Nirb_norect, alpha=0.2)
    ax[0].plot(Ns, mins2Nirb_rect, c="tab:orange", label=r"$\Vert u_h - u_{hH}^N\Vert_{L^2}$ w/ rect.")
    ax[0].plot(Ns, maxs2Nirb_rect, c="tab:orange")
    ax[0].plot(Ns, moy2Nirb_rect, c="tab:orange")
    ax[0].fill_between(Ns, mins2Nirb_rect, maxs2Nirb_rect, alpha=0.2)
    ax[0].axhline(mins2Int_norect[0], color="tab:green", linestyle="--", label=r"$\Vert u_h - u_H\Vert_{L^2}$")
    ax[0].axhline(maxs2Int_norect[0], color="tab:green", linestyle="--")
    ax[0].axhline(moy2Int_norect[0], color="tab:green", linestyle="--")
    ax[0].set_title(r"Error $L^2$")
    ax[0].set_xlabel("N")
    ax[0].set_ylabel("Error")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(Ns[0], Ns[-1])
    ax[0].legend()

    ax[1].plot(Ns, minsInfNirb_norect, c="tab:blue", label=r"$\Vert u_h - u_{hH}^N\Vert_{\infty}$ w/o rect.")
    ax[1].plot(Ns, maxsInfNirb_norect, c="tab:blue")
    ax[1].plot(Ns, moyInfNirb_norect, c="tab:blue")
    ax[1].fill_between(Ns, minsInfNirb_norect, maxsInfNirb_norect, alpha=0.2)
    ax[1].plot(Ns, minsInfNirb_rect, c="tab:orange", label=r"$\Vert u_h - u_{hH}^N\Vert_{\infty}$ w/ rect.")
    ax[1].plot(Ns, maxsInfNirb_rect, c="tab:orange")
    ax[1].plot(Ns, moyInfNirb_rect, c="tab:orange")
    ax[1].fill_between(Ns, minsInfNirb_rect, maxsInfNirb_rect, alpha=0.2)
    ax[1].axhline(minsInfInt_norect[0], color="tab:green", linestyle="--", label=r"$\Vert u_h - u_{hH}^N\Vert_{\infty}$")
    ax[1].axhline(maxsInfInt_norect[0], color="tab:green", linestyle="--")
    ax[1].axhline(moyInfInt_norect[0], color="tab:green", linestyle="--")
    ax[1].set_title(r"Error $L^\infty$")
    ax[1].set_xlabel("N")
    ax[1].set_ylabel("Error")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(Ns[0], Ns[-1])
    ax[1].legend()

    # Save
    tikzplotlib.save("plot.tex")
    plt.show()

# plot_error('/data/scratch/saigre/feel-mbda/nirb/heat/np_1', [2, 3, 4, 5, 10, 15, 20, 25, 50, 100, 175, 200])
plot_error('/data/home/elarif/feelppdb/nirb/heat/np_1', [1, 2, 4, 6, 10, 12, 14, 16, 20, 25, 30, 35, 40, 45, 50, 70, 80, 100])


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
