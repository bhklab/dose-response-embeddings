import matplotlib.pyplot as plt
import seaborn as sns
import umap
import pandas as pd
from damply import dirs
import os 
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import re
root = dirs.PROCDATA / "NCI60" 

def process_atc_code(code_string: str, level: int = 4) -> str:
    """Return the requested ATC hierarchy level (for example N03AX09 -> N03A)."""
    result = re.findall(r"[a-zA-Z]|[^a-zA-Z]+", code_string)
    return "".join(result[:level])


colData = pd.read_csv(dirs.RAWDATA / "HDD/colData.tsv",
                      sep="\t",
                      usecols =  ["HDD.Compound.ID", "Pubchem.CID", "Mechanism.of.Action","ATC.Code"])
files = [
    f for f in os.listdir(root) 
    if os.path.isfile(os.path.join(root, f))
]


reducer = umap.UMAP(random_state=42)

for file in files:
    if file.startswith("."):
        continue
    cancer_type = file.split("_")[0]
    # print(file)
    embeddings = pd.read_csv(dirs.PROCDATA/"NCI60"/file,index_col=0)
    embeddings = embeddings.merge(colData,left_on="PubChem.CID", right_on="Pubchem.CID")
    embeddings = embeddings.dropna(subset=["Mechanism.of.Action","ATC.Code"])
    embeddings['ATC.Code']=[process_atc_code(x,1) for x in embeddings['ATC.Code']]
    X = embeddings[["HS", "E0",  "E_inf", "EC50"]].values 
    X = StandardScaler().fit_transform(X)
    UM = reducer.fit_transform(X)
    plot_df = pd.DataFrame({'UMAP-1':UM[:,0], 
                        'UMAP-2':UM[:,1],
                        'ATC':embeddings['ATC.Code'].values})
    ax =  sns.scatterplot(plot_df,x='UMAP-1',y='UMAP-2',hue='ATC',legend=True)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))
    plt.title(f"{cancer_type} Embeddings vs. ATC Code")
    plt.tight_layout()
    plt.savefig(dirs.RESULTS / f"{cancer_type}_vs_ATC.png")
    plt.close()

    plot_df = pd.DataFrame({'UMAP-1':UM[:,0], 
                            'UMAP-2':UM[:,1],
                            'MOA':embeddings['Mechanism.of.Action'].values})
    ax =  sns.scatterplot(plot_df,x='UMAP-1',y='UMAP-2',hue='MOA',legend=False)
    # sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))
    plt.title(f"{cancer_type} Embeddings vs. Mechanism of Action")
    plt.tight_layout()
    plt.savefig(dirs.RESULTS / f"{cancer_type}_vs_MOA.png")
    plt.close()

    
    