import numpy as np
import pandas as pd
import os.path as path
import re

root = "./data"
chem_gene = pd.read_csv(path.join(root, "raw_data", "CTD_chem_gene_ixns.csv"),
                        header=None,
                        skiprows=range(29),
                        usecols=[0, 3, 8]).to_numpy()

chem_map, gene_map, comp_map, go_map, rela_map = dict(), dict(), dict(), dict(), dict()
adj_chemgo, adj_chemgene, adj_compgo, adj_compgene = set(), set(), set(), set()
adj_chemgo_num, adj_chemgene_num, adj_compgo_num, adj_compgene_num = set(), set(), set(), set()
processed = []
errored = []

for i in range(165):
    try:
        c, g, r = chem_gene[i]

        # ################################################################
        # chemical -- gene process
        # ################################################################
        b = re.match(r"(.*)\[([^\[\]]*)\](.*)", r)  # [chem - relation], [GO]
        if b:
            ind_c, ind_g = chem_map.get(c), go_map.get(b.group(2))
            if ind_c is None:
                chem_map[c] = len(chem_map)
                ind_c = chem_map[c]
            if ind_g is None:
                go_map[b.group(2)] = len(go_map)
                ind_g = go_map[b.group(2)]
            r = re.match(r"{} (.*) ".format(c), b.group(1))
            r = r.group(1).split(' and ')
            for j in r:
                ind_r = rela_map.get(j)
                if ind_r is None:
                    rela_map[j] = len(rela_map)
                    ind_r = rela_map[j]
                adj_chemgo.add((c, b.group(2), j))
                adj_chemgo_num.add((ind_c, ind_g, ind_r))
            processed.append(i)
            continue

        # ################################################################
        # chemical process -- gene
        # ################################################################
        b = re.match(r"\[([^\[\]]*)\](.*)", r)  # [chem proc] - [relation, gene]
        if b:
            ind_c, ind_g = comp_map.get(b.group(1)), gene_map.get(g)
            if ind_c is None:
                comp_map[b.group(1)] = len(comp_map)
                ind_c = comp_map[b.group(1)]
            if ind_g is None:
                gene_map[g] = len(gene_map)
                ind_g = gene_map[g]
            r = re.match(r" (.*) {}".format(g), b.group(2))
            r = r.group(1).split(' and ')
            for j in r:
                ind_r = rela_map.get(j)
                if ind_r is None:
                    rela_map[j] = len(rela_map)
                    ind_r = rela_map[j]
                adj_compgene.add((b.group(1), g, j))
                adj_compgene_num.add((ind_c, ind_g, ind_r))
            processed.append(i)
            continue

        # ################################################################
        # chemical process -- gene process
        # ################################################################
        b = re.match(r"\[([^\[\]]*)\] (.*) \[([^\[\]]*)\]",
                     r)  # [chem - relation], [GO]
        if b:
            ind_c, ind_g = comp_map.get(b.group(1)), go_map.get(b.group(3))
            if ind_c is None:
                comp_map[b.group(1)] = len(comp_map)
                ind_c = comp_map[b.group(1)]
            if ind_g is None:
                go_map[b.group(3)] = len(go_map)
                ind_g = go_map[b.group(3)]
            r = b.group(2).split(' and ')
            for j in r:
                ind_r = rela_map.get(j)
                if ind_r is None:
                    rela_map[j] = len(rela_map)
                    ind_r = rela_map[j]
                adj_compgo.add((b.group(1), b.group(3), j))
                adj_compgo_num.add((ind_c, ind_g, ind_r))
            processed.append(i)
            continue

        # ################################################################
        # chemical -- gene
        # ################################################################
        ind_c, ind_g = chem_map.get(c), gene_map.get(g)
        if ind_c is None:
            chem_map[c] = len(chem_map)
            ind_c = chem_map[c]
        if ind_g is None:
            gene_map[g] = len(gene_map)
            ind_g = gene_map[g]
        r = re.match(r"{} (.*) {}".format(c, g), r)
        r = r.group(1).split(' and ')
        for j in r:
            ind_r = rela_map.get(j)
            if ind_r is None:
                rela_map[j] = len(rela_map)
                ind_r = rela_map[j]
            adj_chemgene.add((c, g, j))
            adj_chemgene_num.add((ind_c, ind_g, ind_r))
        processed.append(i)
    except AttributeError:
        errored.append(i)


chem_set, gene_set, comp_set, go_set, rela_set = set(), set(), set(), set(), set()
adj_chemgo, adj_chemgene, adj_compgo, adj_compgene = set(), set(), set(), set()
processed = []
errored = []

for i in range(165):
    try:
        c, g, r = chem_gene[i]

        # ################################################################
        # chemical --> gene process                 0
        # ################################################################
        b = re.match(r"(.*)\[([^\[\]]*)\](.*)", r)  # [chem - relation], [GO]
        if b:
            r = re.match(r"{} (.*) ".format(c), b.group(1))
            r = r.group(1).split(' and ')
            for j in r:
                adj_chemgo.add((c, b.group(2), j))
            chem_set.add(c)
            go_set.add(b.group(2))
            processed.append(i)
            continue

        # ################################################################
        # chemical process -- gene                     84
        # ################################################################
        b = re.match(r"\[([^\[\]]*)\](.*)", r)  # [chem proc] - [relation, gene]
        if b:
            ind_c, ind_g = comp_map.get(b.group(1)), gene_map.get(g)
            if ind_c is None:
                comp_map[b.group(1)] = len(comp_map)
                ind_c = comp_map[b.group(1)]
            if ind_g is None:
                gene_map[g] = len(gene_map)
                ind_g = gene_map[g]
            r = re.match(r" (.*) {}".format(g), b.group(2))
            r = r.group(1).split(' and ')
            for j in r:
                ind_r = rela_map.get(j)
                if ind_r is None:
                    rela_map[j] = len(rela_map)
                    ind_r = rela_map[j]
                adj_compgene.add((b.group(1), g, j))
                adj_compgene_num.add((ind_c, ind_g, ind_r))
            processed.append(i)
            continue

        # ################################################################
        # chemical process -- gene process
        # ################################################################
        b = re.match(r"\[([^\[\]]*)\] (.*) \[([^\[\]]*)\]",
                     r)  # [chem - relation], [GO]
        if b:
            ind_c, ind_g = comp_map.get(b.group(1)), go_map.get(b.group(3))
            if ind_c is None:
                comp_map[b.group(1)] = len(comp_map)
                ind_c = comp_map[b.group(1)]
            if ind_g is None:
                go_map[b.group(3)] = len(go_map)
                ind_g = go_map[b.group(3)]
            r = b.group(2).split(' and ')
            for j in r:
                ind_r = rela_map.get(j)
                if ind_r is None:
                    rela_map[j] = len(rela_map)
                    ind_r = rela_map[j]
                adj_compgo.add((b.group(1), b.group(3), j))
                adj_compgo_num.add((ind_c, ind_g, ind_r))
            processed.append(i)
            continue

        # ################################################################
        # chemical -- gene
        # ################################################################
        ind_c, ind_g = chem_map.get(c), gene_map.get(g)
        if ind_c is None:
            chem_map[c] = len(chem_map)
            ind_c = chem_map[c]
        if ind_g is None:
            gene_map[g] = len(gene_map)
            ind_g = gene_map[g]
        r = re.match(r"{} (.*) {}".format(c, g), r)
        r = r.group(1).split(' and ')
        for j in r:
            ind_r = rela_map.get(j)
            if ind_r is None:
                rela_map[j] = len(rela_map)
                ind_r = rela_map[j]
            adj_chemgene.add((c, g, j))
            adj_chemgene_num.add((ind_c, ind_g, ind_r))
        processed.append(i)
    except AttributeError:
        errored.append(i)

# ################################################################
# save to file
# ################################################################
# chemical -- gene process
pd.DataFrame(adj_chemgo,
             columns=["ChemicalName", "GeneProcess", "Interaction"]
             ).to_csv(path.join(root, "text_data", "Chemical-GeneProcess.csv"))
pd.DataFrame(adj_chemgo_num,
             columns=["ChemicalIndex", "GeneProcessIndex", "InteractionIndex"]
             ).to_csv(path.join(root, "num_data", "Chemical-GeneProcess.csv"))
# chemical process -- gene
pd.DataFrame(adj_compgene,
             columns=["ChemicalProcess", "GeneSymbol", "Interaction"]
             ).to_csv(path.join(root, "text_data", "ChemicalProcess-Gene.csv"))
pd.DataFrame(adj_compgene_num,
             columns=["ChemicalProcessIndex", "GeneIndex", "InteractionIndex"]
             ).to_csv(path.join(root, "num_data", "ChemicalProcess-Gene.csv"))
# chemical process -- gene process
pd.DataFrame(adj_compgo,
             columns=["ChemicalProcess", "GeneProcess", "Interaction"]
             ).to_csv(path.join(root, "text_data", "ChemicalProcess-GeneProcess.csv"))
pd.DataFrame(adj_compgo_num,
             columns=["ChemicalProcessIndex", "GeneProcessIndex", "InteractionIndex"]
             ).to_csv(path.join(root, "num_data", "ChemicalProcess-GeneProcess.csv"))
# chemical -- gene
pd.DataFrame(adj_chemgene,
             columns=["ChemicalName", "GeneSymbol", "Interaction"]
             ).to_csv(path.join(root, "text_data", "Chemical-Gene.csv"))
pd.DataFrame(adj_chemgene_num,
             columns=["ChemicalIndex", "GeneIndex", "InteractionIndex"]
             ).to_csv(path.join(root, "num_data", "Chemical-Gene.csv"))