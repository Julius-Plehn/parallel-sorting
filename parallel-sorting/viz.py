import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns

# sns.set_theme()
mpl.rcParams["figure.dpi"] = 200


def plot_speedup_2(df):
    ax = df.plot(
        x="# of Processes", y="Speedup", xticks=df["# of Processes"], legend=False
    )
    ax2 = ax.twinx()
    df.plot(
        x="# of Processes",
        y="Time",
        xticks=df["# of Processes"],
        ax=ax2,
        legend=False,
        color="r",
    )
    ax.figure.legend()
    plt.savefig("speedup.pdf")


def plot_speedup(df):
    fig, ax = plt.subplots(figsize=[12, 6])
    sns.pointplot(data=df, x="# of Processes", y="Speedup")
    # plt.legend()
    plt.savefig("speedup.pdf")
    return fig, ax


speedup_df = pd.read_pickle("speedup.pkl")
print(speedup_df)
plot_speedup_2(speedup_df)
