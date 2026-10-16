#!/bin/env python3

import argparse
import os
import sys
import glob

class Gene:
    __slots__ = ("id", "name", "score", "nperm", "nge", "emp", "fdr")
    def __init__(self, id: str, name: str, score: float) -> None:
        self.id    = id
        self.name  = name
        self.score = score
        self.nperm = 0     # permutations in which this gene had a score
        self.nge   = 0     # permutations with a score >= the observed score
        self.emp   = 1.0   # empirical (permutation) p-value
        self.fdr   = 1.0   # Benjamini-Hochberg corrected p-value

perm_path  = ''
obsfile    = ''
outfile    = "fdr_corrected_gene_pvals.tsv"
genes      = dict() # gene_id -> Gene

def get_arguments() -> int:
    """function to get the users input"""

    global perm_path, obsfile, outfile

    # help messages & description
    desc  = "FDR correct each gene's observed P-value using precomputed permutations"
    phelp = "Path to directory containing permutated results"
    Ohelp = "Tsv file containing the observed scores for each gene"
    ohelp = "name of output [optional]"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-p", "--perms", type=str, help=phelp, required=True)
    parser.add_argument("-g", "--genes", type=str, help=Ohelp, required=True)
    parser.add_argument("-o", "--out",   type=str, help=ohelp, default=outfile)

    args      = parser.parse_args()
    perm_path = args.perms
    obsfile   = args.genes
    outfile   = args.out

    assert os.path.isdir(perm_path), f"Unable to find directory {perm_path}"
    assert os.path.isfile(obsfile), f"Unable to find file {obsfile}"

    return 0

def load_genes() -> int:
    """load in the observed genes and their scores"""

    global genes, obsfile

    #
    # expected struct
    #[GeneID][GeneName][TranscriptID][NumSites][AvgScore]
    #

    fh = open(obsfile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields  = line.split('\t')
        assert len(fields) == 5, f"Expected 5 columns in {obsfile}. Offending line:\n{line}"
        gene_id = fields[0]
        name    = fields[1] 
        score   = float(fields[4].strip())
        genes[gene_id] = Gene(gene_id, name, score)

    fh.close()

    if (len(genes) == 0):
        sys.exit(f"Error: Loaded 0 genes from {obsfile}")
    print(f"Loaded {len(genes)} genes")

    return 0

def get_perm_files() -> list[str]:
    """helper function to get the permutation files in the current directory"""

    global perm_path

    perm_files = list()

    if (perm_path[-1] == '/'):
        perm_files = glob.glob(f"{perm_path}*.tsv")
    else:
        perm_files = glob.glob(f"{perm_path}/*.tsv")

    if (len(perm_files) == 0):
        sys.exit(f"Did not find any *.tsv files in {perm_path}")

    print(f"Located {len(perm_files)} permutation files")

    return perm_files

def load_permutations() -> int:
    """add the gene scores into each gene object"""

    global genes

    perm_files = get_perm_files()

    for perm in perm_files:
        fh = open(perm, 'r')

        for line in fh:
            if (len(line) == 0 or line[0] == '#'):
                continue
            fields = line.split('\t')
            assert len(fields) == 5, f"Expected 5 columns in {perm}. Offending line:\n{line}"
            gene_id = fields[0]
            if (gene_id not in genes):
                continue
            score = float(fields[4].strip())
            genes[gene_id].nperm += 1
            if (score >= genes[gene_id].score):
                genes[gene_id].nge += 1

        fh.close()

    return 0

def calc_fdr() -> int:
    """assign a corrected p-value each gene's after multiple testing"""

    global outfile, genes

    # step one, claculate empirical p-value
    tested = list()
    for gene in genes.values():
        if (gene.nperm == 0):
            continue # not in the permutation 
        gene.emp = (gene.nge + 1) / (gene.nperm + 1)
        tested.append(gene)

    if (len(tested) == 0):
        sys.exit(f"Error: No observed genes were found in the permutations")

    # sorted from smallest to largest empirical p-val
    tested.sort(key = lambda gene: gene.emp) 
    m    = len(tested)
    prev = 1.0

    # Benjamini-Hochberg correction
    for i in range(m - 1, -1, -1):
        q             = tested[i].emp * m / (i + 1)
        prev          = min(prev, q)
        tested[i].fdr = prev

    # write the results
    ofh = open(outfile, 'w')
    ofh.write("#GeneID\tGeneName\tObsScore\tNumPerms\tNumAsExtreme\tEmpPval\tFDR\n")

    for gene in tested:
        ofh.write(f"{gene.id}\t{gene.name}\t{gene.score}\t{gene.nperm}\t{gene.nge}\t{gene.emp:.6g}\t{gene.fdr:.6g}\n")
 
    ofh.close()


    ofh.close()

    return 0

def main() -> int:
    """entry point to this application"""

    # get the inputs
    get_arguments()

    # load the observed genes
    load_genes()

    # load the permutations
    load_permutations()

    # now to write each gene's fdr score
    calc_fdr()

    return 0

if __name__ == "__main__":
    main()
