import pandas as pd

debug = True

############################################
# Finding the interactions with seed genes #
############################################

# Read the tables from exercise 2
all_table = pd.read_csv("all.tsv", sep="\t")

# Select he uniprot_id and removing _HUMAN and spaces 
seed_genes = list(map(lambda x: x.strip().replace('_HUMAN', ''), all_table['uniprot_id'].tolist()))

if debug:
    print()
    print("Seed genes:", seed_genes)
    print()

# Read biogrid_human data with interactions
df = pd.read_csv('biogrid_human.txt', sep='\t')

if debug:
    print(df.columns.values)
    print()

# Selecting only the following columns as are the only ones with genes names
ia, ib, ca, cb, aa, ab = 'INTERACTOR_A', 'INTERACTOR_B','OFFICIAL_SYMBOL_A','OFFICIAL_SYMBOL_B','ALIASES_FOR_A','ALIASES_FOR_B'
df2 = df[[ia, ib, ca, cb, aa, ab]]

if debug:
    df5 = df2.loc[(df2[ca].isin(seed_genes))&(df2[cb].isin(seed_genes))]
    print("The number of binary interactions between seed genes is {}".format(df5.shape[0]))
    print()

# Create regular expression to look matches in ALIASES_*
seed_patterns = '|'.join(seed_genes)

# Query to select all rows that have a match with seed_genes names or ALIASES.
df3 = df2.loc[df2[ca].isin(seed_genes)|(df2[cb].isin(seed_genes))|(df2[ia].isin(seed_genes))|(df2[ib].isin(seed_genes))|(df2[aa].str.contains(seed_patterns))|(df2[ab].str.contains(seed_patterns))]

if debug:
    print("We found {} interactions with seed genes".format(df3.shape[0]))
    print()

####################################
# Interactions between interactors #
####################################

# We make a set for seed_genes
seed_genes_set = set(seed_genes)

# We make a set for all the genes in the result
all_seed_genes = set(df3[ca].tolist()).union(set(df3[cb].tolist()))

# We make the set NEW_PROTEINS = ALL_GENES - SEED_GENES
new_proteins = list(all_seed_genes.difference(seed_genes_set))

if debug:
    print("We found {} total genes in the interactions table from which {} are new proteins.".format(len(all_seed_genes), len(new_proteins)))
    print()

# Query to select all rows that match protein A AND protein B in new proteins
df4 = df2.loc[df2[ca].isin(new_proteins)&(df2[cb].isin(new_proteins))]

if debug:
    print("Finally, from {} new proteins we have {} binary interactions".format(len(new_proteins), df4.shape[0]))
    print()

#################
# Output of ex3 #
#################
from bioservices.uniprot import UniProt
u = UniProt(verbose=False)
import time

print(new_proteins)

for x in new_proteins:
    print("######### {} #########".format(x))
    print(u.search("{} AND HUMAN".format(x), frmt="tab", columns="entry name, id, genes"))
    print("##################".format(x))
    print(u.search("{}_HUMAN".format(x), frmt="tab", columns="entry name, id, genes"))
    time.sleep(0.5)
   #print(u.search("{}_HUMAN".format(x), frmt="tab", columns="entry name, id, genes"))
    #res = u.search("({}) AND HUMAN".format("+or+".join(map(lambda x: x+"_HUMAN", new_proteins))), frmt="tab", columns="entry name, id, genes")
print("({}) AND HUMAN".format("+or+".join(map(lambda x: x+"_HUMAN", new_proteins))))
