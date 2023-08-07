import pandas as pd
from urllib.request import urlopen


rule all:
    input:
        "../results/surfaceome_analysis/dmel_adhesion_proteome_flyatlas2.csv"


# Auxiliar functions begin

def excel_stream_to_df(url: str, sheet: int) -> pd.DataFrame:
    "Returns a Pandas DataFrame from an Excel file io stream (urllib)"
    with urlopen(url) as r:
        content = r.read()
    return pd.ExcelFile(content).parse(sheet)


def preprocess_flyatlas2_df(df: pd.DataFrame, proteins: str) -> pd.DataFrame:
    """
    Fix column names, keep only the columns with expression data and
    rows overlapping with the proteins list provided
    """
    # Remove annoying whitespaces from column names and gene names
    df.columns = df.columns.str.strip()
    df["FBgn"] = df["FBgn"].str.strip()
    # Store gene FlyBase ids into indexes
    df.set_index("FBgn", inplace=True)
    # Select tissue expression columns
    dff = df.iloc[:, 8:38]
    # Load protein list
    prots = pd.read_csv(proteins, names=["id"])
    dff = dff[dff.index.isin(prots.id)]
    return dff


def download_flyatlas2_data(url: str, sheet: int, proteins: str, output: str):
    """
    Computes Tau values in a tissue and cell type-specific manner.

    Also calculates the correlation between proteomic data and scRNAseq expression levels
    """

    flyatlas = (excel_stream_to_df(url, sheet)
                .pipe(preprocess_flyatlas2_df, proteins)
                )
    flyatlas.to_csv(output)

# Auxiliar functions end


rule download_flyatlas2_data:
    params:
        url="https://motif.mvls.gla.ac.uk/downloads/FlyAtlas2_gene_data.xlsx",
        proteins="../results/surfaceome_analysis/dmel_adhesion_proteome.csv"
    output:
        outfile="../results/surfaceome_analysis/dmel_adhesion_proteome_flyatlas2.csv"
    run:
        download_flyatlas2_data(params.url, 1, params.proteins, output.outfile)
