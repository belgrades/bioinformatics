import pandas as pd
import io


def search_exact_match(u, gene_name, select_all=False):
    """
    Method that finds giving a gene_name the official name and uniprot ac.
    :param u: The Uniprot object from bioservices
    :param gene_name: The gene name to find.
    :param select_all: If True retrieve a list inside the dict with all genes related to gene_name if exact match fails.
    :return: A dictionary with gene_name and uniprot ac.
    """
    first_query = u.search("{}_HUMAN".format(gene_name), frmt="tab", columns="entry name, id, genes")

    if first_query == "":
        second_query = u.search("{} AND HUMAN".format(gene_name), frmt="tab", columns="entry name, id, genes")
        df = pd.read_csv(io.StringIO(second_query), sep="\t", dtype="object")
        df[['genes_format']] = df[["Gene names"]].astype(str)
        df['genes_format'] = df['genes_format'].apply(lambda x: x.upper())

        # Querying selecting only HUMAN
        df = df.loc[(df["genes_format"].apply(lambda y: gene_name.upper() in y)) & (df["Entry name"].apply(lambda y: "_HUMAN" in y))]

        # Now we have all entries with _HUMAN and gene in genes_format
        # If select_all return all else return first

        for idx, row in df.iterrows():
            name, uniprot = row["Entry name"], row["Entry"]
            if not select_all:
                return {gene_name: {'gene symbol': name.replace("_HUMAN", ""), 'uniprot ac': uniprot}}
    else:
        df = pd.read_csv(io.StringIO(first_query), sep="\t", dtype="object")
        for idx, row in df.iterrows():
            name, uniprot = row['Entry name'], row['Entry']

        return {gene_name: {'gene symbol': name.replace("_HUMAN", ""), 'uniprot ac': uniprot}}
    return {gene_name: {'gene_name': None, 'uniprot ac': None}}


def main():
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
    ia, ib, ca, cb, aa, ab = 'INTERACTOR_A', 'INTERACTOR_B', 'OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B', 'ALIASES_FOR_A', 'ALIASES_FOR_B'
    df2 = df[[ia, ib, ca, cb, aa, ab]]

    if debug:
        df5 = df2.loc[(df2[ca].isin(seed_genes)) & (df2[cb].isin(seed_genes))]
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
    df4 = df2.loc[df2[ca].isin(new_proteins) & (df2[cb].isin(new_proteins))]

    if debug:
        print("Finally, from {} new proteins we have {} binary interactions".format(len(new_proteins), df4.shape[0]))
        print(df4.head(5))
        print()

    #################
    # Output of ex3 #
    #################
    from bioservices.uniprot import UniProt
    u = UniProt(verbose=False)
    import time

    print(new_proteins)
    # Finding
    dictionary_all = {}
    goal = len(new_proteins)
    print("{} proteins to parse".format(goal))
    list_of_nones = []
    for idx, new_gene in enumerate(new_proteins):
        print("\r {}% Processing: {}".format(idx *1.0 / goal * 100.0, new_gene), end="")

        result = search_exact_match(u, new_gene)
        if result[new_gene]['uniprot ac'] is None:
            print("[New None] {}".format(new_gene))
            list_of_nones.append(new_gene)
        else:
            dictionary_all.update(result)
        time.sleep(0.5)
    print(len(dictionary_all), goal)
    print(list_of_nones)


if __name__ == "__main__":
    main()
