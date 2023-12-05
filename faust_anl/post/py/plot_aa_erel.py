import pandas as pd
import seaborn as sns
import time
# import matplotlib.pyplot as plt

if __name__ == "__main__":
    start = time.time()
    input_file_path = "d:tamu_data/exp/si28_c_35/anl/basic_erel.csv"
    mixed_input_file_path = "d:tamu_data/exp/si28_c_35/anl/exactly_two_alphas/basic_erel.csv"
    # plt.rcParams['text.usetex'] = True
    # plt.rcParams['font.family'] = 'serif'
    df = pd.read_csv(input_file_path, chunksize=1000)
    mixed_df = pd.read_csv(mixed_input_file_path, chunksize=1000)
    # df = pd.DataFrame.from_records( pd.read_json(input_file_path)['two_alphas_erel'])

    print(df)

    ax = sns.histplot(df['two_alphas_erel'], bins=10000,
                      element='step', fill=False, color='black')

    sns.histplot(mixed_df['two_alphas_erel'], bins=10000,
                 element='step', fill=False, color='red')

    ax.set_xlim(-.5, 40)
    ax.set(xlabel=r'$2\alpha$ $E_{rel.}$ (MeV)', ylabel=r'Yield')
    fig = ax.get_figure()
    fig.savefig('two_alpha_erel.png')

    stop = time.time()

    print("Plotting took {}".format(stop - start))
