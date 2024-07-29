import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


data = pd.read_csv("../data/quinone_metabolism.csv", index_col=0)
seed = 42
n_permutations = 10000


random_state = np.random.RandomState(seed=seed)


def compute_statistics(data):
    #
    perc = (
        data.loc[data["Q_type_tested"] == "UQ-O2dep"].groupby("metabolism").count() /
        data.loc[data["Q_type_tested"] == "UQ-O2dep"].shape[0])
    uq_O2dep_perc_aero = perc.loc["aerobic"].values[0]

    perc = (
        data.loc[data["Q_type_tested"] == "UQ-O2dep_indep"].groupby("metabolism").count() /
        data.loc[data["Q_type_tested"] == "UQ-O2dep_indep"].shape[0])
    uq_O2dep_indep_perc_aero = perc.loc["aerobic"].values[0]

    perc = (
        data.loc[data["Q_type_tested"] == "HP+LP"].groupby("metabolism").count() /
        data.loc[data["Q_type_tested"] == "HP+LP"].shape[0])
    hplp_perc_aero = perc.loc["aerobic"].values[0]

    perc = (
        data.loc[data["Q_type_tested"] == "LP"].groupby("metabolism").count() /
        data.loc[data["Q_type_tested"] == "LP"].shape[0])
    try:
        lp_perc_aero = perc.loc["aerobic"].values[0]
    except KeyError:
        lp_perc_aero = 0

    return min(
        uq_O2dep_perc_aero - uq_O2dep_indep_perc_aero,
        uq_O2dep_indep_perc_aero - hplp_perc_aero,
        hplp_perc_aero - lp_perc_aero
        )

statistic = compute_statistics(data)

permutated_statistics = np.zeros(n_permutations+1)

for permutation_id in range(n_permutations):
    data.loc[:, "metabolism"] = random_state.permutation(data["metabolism"])
    permutated_statistics[permutation_id] = compute_statistics(data)


# Adding the original statistic to the list of permutation
permutated_statistics[-1] = statistic

pvalue = np.mean(statistic <= permutated_statistics)
print(pvalue)


fig, ax = plt.subplots(figsize=(10, 6), tight_layout=True)
ax.hist(permutated_statistics, bins=50, color="0")
ax.axvline(statistic, linewidth=2)
ax.spines["right"].set_linewidth(0)
ax.spines["top"].set_linewidth(0)
ax.set_title("(pval: %0.2f)" % pvalue, fontweight="bold")
