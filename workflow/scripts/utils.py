import gzip
import requests
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from urllib.request import urlopen


def get_url_links(url: str) -> list:
    """
    Returns all the URL links (markerd with an <a> tag) inside HTMl content
    """

    # Load the HTML file
    with requests.get(url) as r:
        soup = BeautifulSoup(r.content, "html.parser")
    # Extract all hyperlinks
    links = [link.contents[0]
             for link in soup.find_all("a")]
    return links


def download_compressed_file(url: str, outfile: str) -> None:
    with requests.get(url, stream=True) as r:
        # with gzip.GzipFile(fileobj=r) as uncompressed:
        with open(outfile, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8192):
                fh.write(chunk)
    return None


def tau_index(v: pd.Series) -> pd.Series:
    return (np.sum(1 - v/np.max(v)) / (len(v) - 1))


def calculate_tau(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Tau index for the columns matrix/DataFrame
    """

    return df.apply(tau_index, axis=1)


### FlyAtlas2 and scRNAseq expression specificity analysis


def aggregate_variable(df: pd.DataFrame, groups: list, variable: str) -> pd.DataFrame:
    """
    Aggregates DataFrame variable mean value for the groups provided
    """
    
    return df.groupby(groups)[variable].mean().reset_index()


def average_replicates_expression(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns replicates mean from an expression matrix
    """

    cells = ["FBgn", "Subtype", "Replicate", "Timepoint"]
    replicates = ["FBgn", "Subtype", "Timepoint"]
    return (df
            .pipe(aggregate_variable, cells, "Expression")
            .pipe(aggregate_variable, replicates, "Expression")
            )


def select_timepoint(df: pd.DataFrame, timepoint:str) -> pd.DataFrame:
    """
    Return only data beloging to the timepoint provided
    """

    dff = pd.pivot(df[df["Timepoint"] == timepoint],
                   index="FBgn",
                   columns="Subtype",
                   values="Expression"
                   )

    return dff


def cell_specificity(df: pd.DataFrame, timepoint: str) -> pd.DataFrame:
    """
    Tau expression specificy for each cell type at the given timepoint
    """

    return (df
            .pipe(average_replicates_expression)
            .pipe(select_timepoint, timepoint)
            .pipe(calculate_tau)
            )


def merge_t4t5_celltypes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Creates a fictional 'T4T5' celltype
    """

    df["T4T5"] = df.loc[:, ["T4a", "T4b", "T4c", "T4d",
                            "T5a", "T5b", "T5c", "T5d"]].max(axis=1)
    return df


def intersect_scrnaseq_proteome(df, proteome):
    """
    Keep only scRNAseq-expressed genes that were also detected on the proteomic dataset
    """

    prt = pd.read_csv(proteome, index_col="Protein")
    prt = prt[(prt.T4_signal > 0.18) & (prt.T5_signal > 0.18)].mean(axis=1)
    dff = pd.concat([df, prt], axis=1).dropna()
    
    return dff


def scrnaseq_proteome_correlation(df, proteome):
    """
    """

    dff = (df
           .pipe(average_replicates_expression)
           .pipe(intersect_scrnaseq_proteome, proteome)
            )
    
    for cell in ["T4a", "T4b", "T4c", "T4d", "T5a", "T5b", "T5c", "T5d"]:
        prt = pd.read_csv(proteome, index_col="Protein")
        prt = prt[(prt.T4_signal > 0.18) & (prt.T5_signal > 0.18)].mean(axis=1)
        dff2 = dff[dff["Subtype"] == cell]
        dff2 = dff2[dff2.index.isin(prt.index)]
        print(dff2.groupby(["FBgn"])["Expression"].std())
        # sns.boxplot(x="FBgn", y="Expression", hue="Timepoint", data=dff2)
        # plt.savefig(f"{cell}_expression.svg")
        # plt.close()
    

    for timepoint in ["24h", "36h", "48h", "60h", "72h", "84h", "96h"]:
        rows = []
        tmp = (dff
               .pipe(select_timepoint, timepoint)
               # .pipe(merge_t4t5_celltypes)
               .pipe(intersect_scrnaseq_proteome, proteome)
               )
        for col in tmp:
            rows.append([str(col), np.mean(tmp[col])])
            # print(f"{timepoint}, {col}, {np.mean(tmp[col])}, {spearmanr(tmp[col], tmp[0]).correlation}")
        ints = pd.DataFrame(rows, columns=["Celltype", "Intensity"]).sort_values("Intensity", ascending=False)
        print(ints.iloc[:10, :])
        sns.scatterplot(y=ints["Intensity"], x=ints["Celltype"])
        plt.savefig(f"{timepoint}_intensities.svg")
        plt.close()
    return dff


def excel_stream_to_df(url: str, sheet: str) -> pd.DataFrame:
    """
    Returns a Pandas DataFrame from an Excel file io stream (urllib)
    """

    with urlopen(url) as r:
        content = r.read()
    return pd.ExcelFile(content).parse(sheet)


def preprocess_flyatlas2_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fix column names and filter only the expression columns from FlyAtlas2 data
    """

    # Remove annoying whitespaces from column names and gene names
    df.columns = df.columns.str.strip()
    df["FBgn"] = df["FBgn"].str.strip()
    # Store gene FlyBase ids into indexes
    df.set_index("FBgn", inplace=True)
    # Select tissue expression columns and compute Tau
    return df.iloc[:, 8:38]


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
    plt.savefig("tissue_specificity_plot.svg")

    return p


def expression_specificity_analysis(url: str, sheet: str, scrnaseq: str, proteome: str):
    """
    Computes Tau values in a tissue and cell type-specific manner.

    Also calculates the correlation between proteomic data and scRNAseq expression levels
    """

    # tau_tissue = (excel_stream_to_df(url, sheet)
    #               .pipe(preprocess_flyatlas2_df)
    #               .pipe(calculate_tau)
    #               )

    # tau_cell = (pd.read_csv(scrnaseq)
    #             .pipe(cell_specificity, "96h")
    #             )

    cell_exp = (pd.read_csv(scrnaseq)
                .pipe(scrnaseq_proteome_correlation, proteome)
                )
