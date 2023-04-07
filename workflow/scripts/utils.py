import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from urllib.request import urlopen


def urllib_request(url: str) -> bytes:
    """
    Returns the content from an URL link
    """


    with urlopen(url) as r:
        content = r.read()

    return content


### FlyAtlas2 data analysis utilities

def excel_stream_to_df(url: str, sheet: str) -> pd.DataFrame:
    """
    Returns a Pandas DataFrame from an Excel file io stream (urllib)
    """

    return pd.ExcelFile(urllib_request(url)).parse(sheet)


def tau_index(v: pd.Series) -> pd.Series:
    """
    Returns Tau expression index from a vector/series
    """

    return (np.sum(1 - v/np.max(v)) / (len(v) - 1))


def preprocess_flyatlas2_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fix column names and filter only the expression columns from FlyAtlas2 data
    """

    # Remove annoying whitespaces from columns
    df.columns = df.columns.str.strip()
    # Store gene FlyBase ids into indexes
    df.set_index("FBgn", inplace=True)
    # Select tissue expression columns and compute Tau
    return df.iloc[:, 8:38]


def calculate_tau(df: pd.DataFrame) -> pd.Series:
    """
    Compute Tau index for a matrix/DataFrame
    """

    return df.apply(tau_index, axis=1)


def plot_flyatlas2_tau(data: pd.Series) -> plt.Axes:
    """
    Plots an histogram od Tau (tissue-specific expression) distribution
    """

    # Figure general parameters
    sns.set(style="darkgrid",
            font_scale=1.8,
            rc={"figure.figsize":(9, 9),
                "lines.linewidth": 3})
    # Plot data
    p = sns.histplot(data,
                     color="black")
    # Set labels legand
    p.set_xlabel("Tau")
    p.set_ylabel("Gene Count")
    # Save figure to file
    plt.tight_layout(h_pad=2)
    plt.savefig("flyatlas2_tau_distribution.svg")

    return p


def flyatlas2_tau_analysis(url: str, sheet: str) -> plt.Axes:
    """
    Runs Tau analysis of FlyAtlas2 data in one step
    """

    df = excel_stream_to_df(url, sheet)

    return (df
            .pipe(preprocess_flyatlas2_df)
            .pipe(calculate_tau)
            .pipe(plot_flyatlas2_tau)
            )
