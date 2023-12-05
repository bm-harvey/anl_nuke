import pandas as pd
import seaborn as sns
import time
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    start = time.time()
    real_file_path = "d:tamu_data/exp/si28_c_35/anl/dd/real/erel.csv"
    mixed_file_path = "d:tamu_data/exp/si28_c_35/anl/dd/mixed/erel.csv"
    real_df = pd.read_csv(real_file_path)
    mixed_df = pd.read_csv(mixed_file_path)

    bins = np.linspace(-5, 50, 300)

    ax = sns.histplot(real_df['e_rel'], bins=bins,
                      element='step', fill=False, color='black', log_scale=(False, False), stat='density')
    ax.set(xlabel=r'dd $E_{rel.}$ (MeV)', ylabel=r'Yield')
    fig = ax.get_figure()
    fig.savefig('dd_erel.png')

    # plt.cla()
    ax = sns.histplot(mixed_df['e_rel'], bins=bins,
                      element='step', fill=False, color='red', stat='density')
    # ax.set_xlim(-.5, 40)
    ax.set(xlabel=r'$2\alpha$ $E_{rel.}$ (MeV)', ylabel=r'Yield')
    fig = ax.get_figure()
    fig.savefig('mixed_dd_erel.png')

    stop = time.time()

    print("Plotting took {}".format(stop - start))
