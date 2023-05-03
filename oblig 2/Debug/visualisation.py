# visualisation.py


# 1: plot correlation vs r(site)
# 2: plot real part of average magnetisation per site <m> vs T in 2D
# TODO 3: plot <|m|^2> vs T for L=16
# TODO 4: plot gamma = m4/m2^2 


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

L = 16

# part 1: correlation
corr_mc_025 = pd.read_csv("correlation,3,16,0.25.csv", header=None, names=["r", "corr"])
corr_mc_05 = pd.read_csv("correlation,3,16,0.50.csv", header=None, names=["r", "corr"])

def corr_exact(r):
    return (eig2**r * eig1**(L-r) + eig1**r * eig2** (L-r) + eig1 ** L )/Z


r = corr_mc_025['r']
plt.scatter(r, corr_mc_025['corr'], label="mc_0.25")
plt.scatter(r, corr_mc_05['corr'], label="mc_0.5")
for T in (0.25, 0.5):
    ebj = np.exp(1/T)
    eig1 = ebj - 1
    eig2 = ebj + 2
    Z = 2* eig1 ** L + eig2 ** L
    plt.plot(r, corr_exact(r), label="exact" + str(T))

    plt.title("Correlations C(r) vs r, L=16, T=" + str(T))
plt.legend()
plt.ylabel("C(r)")
# plt.savefig("../plots/correlation-mc-exact-cmpr.pdf")
plt.close()


# part 2 
for L in (8,16,32):
    mag = pd.read_csv(str(L) + " magnetisations_vs_T.csv")
    if L == 16:
        mag.plot(x="T", y=["m","m2"])#, "m4"])
        plt.title("Magnetisation per Site vs T/J")
        plt.savefig("../plots/average_m.pdf")
        plt.close()

# part 3: the gamma ratio for different L's
    gamma = np.array(mag["m4"] / mag["m2"]**2)
    # print(gamma)
    plt.plot(mag["T"], gamma, label="L=" + str(L))
    plt.ylabel(r"$\gamma$")
    plt.xlabel("T")
    plt.legend()

plt.savefig("../plots/gamma.pdf")
plt.show()