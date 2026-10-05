#!/bin/env python3

import argparse
import os
import gzip
import sys
import matplotlib.pyplot as plt

class Chrom:
    def __init__(self, name: str, length: int) -> None:
        self.name           = name
        self.length         = length
        self.gene_ids       = list()
        self.gene_colors    = list()
        self.gene_values    = list()
        self.gene_positions = list()

    def add_xposition(self, x_pos) -> None:
        for i in range(len(self.gene_positions)):
            self.gene_positions[i] += x_pos

fai_file    = ''
ann_file    = ''
genes_file  = ''
sig_file    = '' # significant
chroms      = dict() # str -> Chrom
chroms_list = list() # tuples (chrom, int) # same as map
min_len     = 0
xpos_map    = dict() # x-axis position to begin plotting
genomeSize  = 0
ymax        = 0

def get_arguments() -> int:
    """function to get the users input"""

    global fai_file, ann_file, genes_file, min_len, sig_file

    # help messages & description
    desc  = "Plot the averaged dxy score per gene across the genome"
    xhelp = "Avg_gene_scores.tsv generated from calc_gene_avg_pixy"
    ahelp = "Annotation file containing the location of the genes"
    fhelp = "fasta index file with the lengths of the chromosomes"
    chelp = "Optional file of gene_ids to color red to distinguish them"
    mhelp = "Minimum length of chromosome to include [default: 0]"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-g", "--genes", required=True,  type=str, help=xhelp)
    parser.add_argument("-a", "--ann",   required=True,  type=str, help=ahelp)
    parser.add_argument("-f", "--fai",   required=True,  type=str, help=fhelp)
    parser.add_argument("-c", "--color", required=False, type=str, help=chelp)
    parser.add_argument("-m", "--min",   type=int, help=mhelp, default=min_len)

    args = parser.parse_args()

    assert os.path.isfile(args.genes),f"Could not locate {args.genes}"
    assert os.path.isfile(args.ann),  f"Could not locate {args.ann}"
    assert os.path.isfile(args.fai),  f"Could not locate {args.fai}"
    assert args.min >= 0, "--min cannot be less than 0"

    genes_file = args.genes
    ann_file   = args.ann
    fai_file   = args.fai
    min_len    = args.min

    if (args.color != None):
        assert os.path.isfile(args.color),  f"Could not locate {args.color}"
        sig_file = args.color

    return 0

def load_chrom_lengths() -> int:
    """read in the chromosome lengths from the index file"""

    global fai_file, chroms, chroms_list, min_len, genomeSize

    fh = open(fai_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        chrom  = fields[0]
        length = int(fields[1])
        if (length >= min_len or chrom == "MT"):
            chroms[chrom] = Chrom(chrom, length)
            genomeSize   += length
            chroms_list.append((chrom, length))

    fh.close()

    if (len(chroms) == 0):
        msg = f"No chromosomes found with a minumum length of {min_len}"
        sys.exit(msg)

    # now set the x positions
    chroms_list.sort(key=lambda x: x[1], reverse=True) # largest chroms first
    xstart = 0
    for chrom, chrom_len in chroms_list:
        xpos_map[chrom] = xstart
        xstart         += chrom_len

    return 0

def load_genes() -> dict[str, float]:
    """parse the tsv file containing the dxy values"""

    global genes_file

    fh        = gzip.open(genes_file, "rt") if genes_file.endswith(".gz") else open(genes_file, 'r')
    genes_map = dict()

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        parts = line.split('\t')
        assert len(parts) == 5, f"Malformed line found in {genes_file}: \n{line}"
        gene_id = parts[0]
        value   = float(parts[4].strip())
        genes_map[gene_id] = value

    fh.close()

    if (len(genes_map) == 0):
        msg = f"No genes found in {genes_file}"
        sys.exit(msg)

    print(f"Loaded {len(genes_map)} genes from {genes_file}")

    return genes_map

def load_colors() -> set[str]:
    """return a set of the gene ids to color red"""

    global sig_file

    if (sig_file == ''):
        return set() # nothing to do hear

    gene_ids = set()
    fh       = gzip.open(sig_file, "rt") if sig_file.endswith(".gz") else open(sig_file, 'r')

    for line in fh:
        if (len(line) == 0):
            continue
        gene_ids.add(line.strip())

    fh.close()

    return gene_ids

def get_gene_id(attrb: str) -> str:

    fields  = attrb.split(';')
    gene_id = ''
    id_     = ''

    if (len(fields) == 1):
        if ("gene_id" in fields[0]):
            gene_id = fields[0].replace("gene_id", '')
            gene_id = gene_id.strip(' "\n')
        elif ("ID" in fields[0]):
            gene_id = fields[0].replace("ID", '')
            gene_id = gene_id.strip(' "\n=')
        else:
            gene_id = fields[0].strip(' "\n')
        return gene_id
    
    for field in fields:
        field = field.strip(' "')
        if ((field.startswith("gene_id") == False) and (field.startswith("ID") == False)):
            continue
        subfields = field.strip(' "\n,=').split(' ')
        if (len(subfields) == 1 and '=' in subfields[0]):
            subfields = field.strip(' "\n,=').split('=')
        record_id = subfields[-1]
        record_id = record_id.strip(' "\n')

        # hold onto this id if no gene_id found
        if (field[0] == 'I'):
            id_     = record_id
        else:
            gene_id = record_id
            
    if (gene_id == ''):
        gene_id = id_ # assume an ID= was found

    return gene_id

def get_gene_coordinates() -> int:
    """load the gene coordinates for the gene_ids in the gene_map"""

    global ann_file, chroms

    # get the gene_ids
    gene_map    = load_genes()
    found_genes = list()
    fh          = gzip.open(ann_file, "rt") if ann_file.endswith(".gz") else open(ann_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 9, f"Malformed line found in {ann_file}: \n{line}"
        if (fields[2] != "gene"):
            continue
        gene_id = get_gene_id(fields[8])
        if (gene_id in gene_map):
            midpoint = (int(fields[3]) + int(fields[4])) / 2
            found_genes.append((fields[0], gene_id, midpoint, gene_map[gene_id]))

    fh.close()

    if (len(found_genes) != len(gene_map)):
        msg = f"Unable to find all genes in {ann_file}"
        sys.exit(msg)
    else:
        print("Found all genes in the annotation file")

    # now to place them onto the chromosomes
    found_genes.sort(key = lambda e: (e[0], e[2]))

    global ymax

    significant = load_colors()

    for entry in found_genes:
        chrom   = entry[0]
        gene_id = entry[1]
        bp_pos  = entry[2]
        value   = entry[3]
        color   = "#E41A1C" if gene_id in significant else "#6FA8DC"
        if (chrom in chroms): # due to min length filter
            chroms[chrom].gene_ids.append(gene_id)
            chroms[chrom].gene_values.append(value)
            chroms[chrom].gene_positions.append(bp_pos)
            chroms[chrom].gene_colors.append(color)
            ymax = max(ymax, value)

    return 0

def plot_manhattan() -> int:
    """Plot all the genes as a scatter plot"""

    global genomeSize, chroms, chroms_list, ymax

    # create a shared drawing surface
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (20, 5))

    ax.set_ylim(0, ymax * 1.05)
    ax.set_xlim(0, genomeSize)

    current_pos = 0
    xlabels     = list()
    xticks      = list()

    for i, chrom_tuple in enumerate(chroms_list, start=1):
        chrom = chroms[chrom_tuple[0]]
        # draw the x-axis chrom label
        xpos = current_pos + (chrom.length / 2)
        xticks.append(xpos)
        xlabels.append(chrom.name)
        chrom.add_xposition(current_pos)
        current_pos += chrom.length
        if (i != len(chroms_list)): # not the last chrom
            ax.axvline(current_pos, color = "black", linestyle = "--")

        # now add the scatter point
        plt.scatter(chrom.gene_positions, chrom.gene_values, color = chrom.gene_colors, s = 5)

    # add chrom labels
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    # add y label
    ax.set_ylabel(r"Average Gene $D_{XY}$")

    plt.savefig("avg_gene_dxy.png", dpi = 500, bbox_inches = "tight")

    return 0

def main() -> int:
    """entry point to this application"""

    # get the input files
    get_arguments()

    # load the chromosome lengths & create chrom objects
    load_chrom_lengths()

    # get the gets
    get_gene_coordinates()

    # now plot
    plot_manhattan()
    
    return 0

if __name__ == "__main__":
    main()
