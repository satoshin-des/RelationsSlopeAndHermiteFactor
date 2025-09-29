import sato_exp
import matplotlib.pyplot as plt
import pandas as pd

N = 100
for N in range(100, 151, 10):
    C = sato_exp.svp_challenge(N, 0)
    print(C)

    print(C.bkz(beta=40, max_loops=10, pruning=True))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("number of tours")
    g = pd.read_csv("hf_log.csv")
    ax.plot([i for i in range(len(g['val']))], g['val'], "--", marker = "", label=r"$\frac{2}{n-1}\log \mathrm{hf}(\boldsymbol{B})$")
    g = pd.read_csv("sl_log.csv")
    ax.plot([i for i in range(len(g['val']))], g['val'], marker = "", label="GSA-slope")

    plt.tick_params()
    plt.legend()
    plt.savefig(f"{N}.png")
