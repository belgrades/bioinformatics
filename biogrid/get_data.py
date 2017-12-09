from bioservices import PSICQUIC
from bioservices import UniProt
import pandas as pd

url_file = "seed_genes.csv"

df_sd = pd.read_csv(url_file)
all_table = pd.read_csv("all.tsv", sep="\t")

s = PSICQUIC(verbose=False)
u = UniProt()

for gene_s in all_table['uniprot_id']:
    gene = gene_s.strip()
    res = u.mapping(fr="ACC", to="BIOGRID_ID", query=gene)
    if(len(res)>0):
        print("Gene", gene)
        print(res)
        print(res[gene])
        print(res[gene][0])
        result = s.query('biogrid', res[gene][0])
        for x in result:
            print("Nuovo")
            for y in x:
                print(y)
        break
    else:
        print(gene)

    print()
