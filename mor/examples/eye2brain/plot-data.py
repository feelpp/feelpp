# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import json
import tikzplotlib

def plot(path):
    print(path)

    with open(path, 'r') as f:
        data = json.load(f)

    N = data["N"]
    algo = data["algo"]
    sampling_size = data["sampling-size"]
    names = data["Names"]

    xs = np.array(list(range(N)))

    firstOrderIndices = data["FirstOrder"]["values"]
    firstOrderIntervals = data["FirstOrder"].get("intervals", None)

    totalOrderIndices = data["TotalOrder"]["values"]
    totalOrderIntervals = data["TotalOrder"].get("intervals", None)


    firstOrderIndices = np.array(firstOrderIndices)
    totalOrderIndices = np.array(totalOrderIndices)
    print("First order sum", firstOrderIndices.sum())
    print("Total order sum", totalOrderIndices.sum())

    if firstOrderIntervals is not None:
        firstOrderIntervals = np.array(firstOrderIntervals)
        totalOrderIntervals = np.array(totalOrderIntervals)

        print("Max error first order", np.max(firstOrderIntervals[:, 1] - firstOrderIntervals[:, 0]))

        plt.errorbar(xs, firstOrderIndices, yerr=[firstOrderIndices - firstOrderIntervals[:,0], firstOrderIntervals[:,1] - firstOrderIndices], fmt='o', label="First order")
        plt.errorbar(xs+0.2, totalOrderIndices, yerr=[totalOrderIndices - totalOrderIntervals[:,0], totalOrderIntervals[:,1] - totalOrderIndices], fmt='o', label="Total order")

    else:
        plt.scatter(xs, firstOrderIndices, label="First order")
        plt.scatter(xs+0.2, totalOrderIndices, label="Total order")

    plt.xticks(xs+0.1, names)

    plt.axhline(y=0, color='r', linestyle='-')
    plt.axhline(y=1, color='r', linestyle='-')

    plt.title(f"Samping of size {sampling_size} ({algo})")
    plt.legend()
    plt.grid()
    tikzplotlib.save("plot.tikz")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) > 1 and 'ip' not in sys.argv[1]:
        path = sys.argv[1]
    else:
        # path = "/data/scratch/saigre/feel-nirb/feelpp_mor_sensitivity_analysis/np_1/B1.json"
        path = "feelpp_mor_sensitivity_analysis/np_1/sensitivity-bootstrap.json"
    plot(path)

# %%
