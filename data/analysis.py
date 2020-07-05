import pandas as pd
import os.path as path

root = "./data/raw_data"
chem = pd.read_csv(path.join(root, "CTD_chemicals.csv"),
                   header=None,
                   skiprows=range(29),
                   usecols=[1])
gene = pd.read_csv(path.join(root, "CTD_genes.csv"),
                   header=None,
                   skiprows=range(29),
                   usecols=[2])
dise = pd.read_csv(path.join(root, "CTD_diseases.csv"),
                   header=None,
                   skiprows=range(29),
                   usecols=[1])
chem_gene = pd.read_csv(path.join(root, "CTD_chem_gene_ixns.csv"),
                        header=None,
                        skiprows=range(29),
                        usecols=[1, 4, 9]).to_numpy()
chem_dise = pd.read_csv(path.join(root, "CTD_chemicals_diseases.csv"),
                        header=None,
                        skiprows=range(29),
                        usecols=[1, 4, 7]).to_numpy()
gene_dise = pd.read_csv(path.join(root, "CTD_genes_diseases.csv"),
                        header=None,
                        skiprows=range(29),
                        usecols=[1, 3, 6]).to_numpy()



