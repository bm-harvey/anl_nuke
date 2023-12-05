import pandas as pd
import seaborn as sns
import time
import matplotlib.pyplot as plt

if __name__ == "__main__":
    start = time.time()
    input_file_path = "../time_to_generate.csv"
    df = pd.read_csv(input_file_path)
    # Using reset_index, inplace=True
    df.reset_index(inplace=True)
    print(df)

    # ax = sns.histplot(df, y='time', x='index', bins=500)
    ax = sns.scatterplot(df, y='time', x='index', color='black', size=.5)

    ax.set(xlabel='Event', ylabel=r'Time to generate event (ns)', yscale="log")
    fig = ax.get_figure()

    fig.savefig('time_to_generate_mixed.png')

    stop = time.time()

    print("Plotting took {}".format(stop - start))
