#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

region  = ''
vcf     = ''
popFile = ''     # 4 column file of [sample][color][pop][sex]
minVar  = 2      # minimum variants
samples = dict() # sample  -> Sample
annPCA  = False  # annotate the PCA plot with the sample's id

class Sample:

    __slots__ = ("id", "color", "population", "sex", "variants")

    def __init__(self, id: str, color: str, population: str, sex: str) -> None:
        self.id         = id
        self.color      = color
        self.population = population
        self.sex        = sex
        self.variants   = list()
    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""

    global minVar, region, vcf, popFile, annPCA
    
    desc  = "Generate a local PCA using variants within a specific region"
    shelp = "Four column tsv file containing sample and hex value for PCA color point, population ID, and sex [M/F]"
    vhelp = "Vcf file containing samples found in --samples"
    rhelp = "region to plot. Format expected CHR:Start-END"
    mhelp = "Minimum number of variant sites required"
    lhelp = "If set, label each point in the PCA with sample IDs"
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-s", "--samples", required=True,  type=str, help=shelp)
    parser.add_argument("-v", "--vcf",     required=True,  type=str, help=vhelp)
    parser.add_argument("-r", "--region",  required=True, type=str, help=rhelp)
    parser.add_argument("-m", "--min",     default=minVar, type=int, help=mhelp)
    parser.add_argument("-l", "--label",   action="store_true", help=lhelp)

    args = parser.parse_args()
    assert os.path.isfile(args.samples), f"Could not find {args.samples}"
    assert os.path.isfile(args.vcf),     f"Could not find {args.vcf}"
    assert ':' in args.region,           f"Region expected format: CHR:Start-END"
    assert args.min > 1,                 f"--min must be at least 2"

    popFile = args.samples
    vcf     = args.vcf
    region  = args.region
    minVar  = args.min
    annPCA  = args.label

    return 0


def load_samples() -> int:
    """load each sample and their color in the pca plots"""

    global samples, popFile

    fh = gzip.open(popFile, "rt") if popFile.endswith(".gz") else open(popFile, 'r')

    for line in fh:
        if (len(line) == 0):
            continue
        fields = line.split('\t')
        assert len(fields) == 4, f"Malformed line detected in {popFile}:\n{line}"
        sample_id  = fields[0]
        color      = fields[1]
        population = fields[2]
        sex        = fields[3].strip()
        assert color[0] == '#', f"Second column in {popFile} should be a hex value. Found {color}"
        assert sex in "MF", "Valid sex info options: M, F"
        samples[sample_id] = Sample(sample_id, color, population, sex)

    fh.close()

    if (len(samples) == 0):
        sys.exit(f"No samples found in {popFile}")

    print(f"Loaded {len(samples)} samples from {popFile}")

    return 0

def parse_region() -> tuple[str, int, int]:
    """parse the region and return the coordinates"""

    global region

    fields = region.split(':')

    assert len(fields) == 2, "Error: Expected format for region CHR:Start-END"

    chrom     = fields[0]
    positions = fields[1].split('-')

    assert len(positions) == 2, "Error: Expected format for region CHR:Start-END"

    left  = int(positions[0])
    right = int(positions[1])

    if (left < 0 or right < 0):
        sys.exit("Negative values for region are not valid")
    if (right < left):
        sys.exit("Start bp cannot be lower than End bp")

    return (chrom, left, right)

def parse_vcf() -> int:
    """go through the vcf and find exonic variants for the target genes"""

    global vcf, samples

    samples_ls = list()
    target_reg = parse_region()
    target_chr = target_reg[0] + '\t'
    left       = target_reg[1]
    right      = target_reg[2]
    sites       = 0
    fh          = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')
    genomap     =  {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2,
                    "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2}

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields     = line.split('\t')
                fields[-1] = fields[-1].strip() # remove end line char
                for i in range(9, len(fields)):
                    sample = fields[i]
                    if (sample not in samples):
                        sys.exit(f"Vcf must have the same samples as --samples. Encountered {sample}")
                    samples_ls.append(sample)
                assert len(samples_ls) == len(samples), f"Not every sample in --samples was found in vcf file"
            continue
        if (len(samples_ls) == 0):
            sys.exit(f"ERROR: Did not encounter #CHROM header in {vcf}")
        elif (line.startswith(target_chr) == False):
            continue # skip other chroms

        head = line.split('\t', 2)
        pos  = int(head[1])
        
        if ((left <= pos <= right) == False):
            continue # different part of the target chromosome

        sites  += 1
        # add each sample's genotype to each gene
        fields     = line.split('\t')
        fields[-1] = fields[-1].strip()
        missing    = list()
        total      = 0
        count      = 0
        for i in range(9, len(fields)):
            sample = samples_ls[i - 9]
            geno   = fields[i].split(':')[0]
            if (geno.count('.') == 0):
                count += 1
                score  = genomap[geno]
                total += score
                samples[sample].variants.append(score)
            else:
                missing.append(sample)
        # treat missing value as averages
        score = 0 if total == 0 else total / count
        for sample in missing:
            samples[sample].variants.append(score)

    fh.close()

    if (sites == 0):
        sys.exit("No variant sites found in any target genes")
    elif (sites < minVar):
        sys.exit(f"{sites} vcf records detected. Minimum required: {minVar}. Exiting..")
    else:
        print(f"Total number of sites found: {sites}")

    return 0

def plot_pca(samples_ls: list[str], scores: np.ndarray, var: np.ndarray, nsites: int) -> int:
    """scatter plot of PC1 and PC2 colored by the user's sample colors"""

    global samples, region, annPCA

    # create the colors and legends features
    colors = list()
    shapes = list()
    legmap = dict()

    for sample in samples_ls:
        color = samples[sample].color
        shape = 'o' if samples[sample].sex == 'M' else '^'
        colors.append(color)
        shapes.append(shape)
        if (color not in legmap):
            legmap[color] = samples[sample].population

    handles = list()
    # create the legend
    for color, legID in legmap.items():
        legend = Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color,
                        markeredgecolor = "black", markeredgewidth = 0.5, markersize = 8, label = legID)
        handles.append(legend)

    # add the sex shapes
    sexes = {'o': "Male", '^': "Female"}
    for shape, label in sexes.items():
        if (shape not in shapes):
            continue
        legend = Line2D([0], [0], marker = shape, color = 'w', markerfacecolor = "gray",
                        markeredgecolor = "black", markeredgewidth = 0.5, markersize = 8, label = label)
        handles.append(legend)

    # now plot
    fig, ax   = plt.subplots(figsize = (6, 6))
    colors_np = np.array(colors)
    shapes_np = np.array(shapes)
    for shape in ('o', '^'):
        mask = shapes_np == shape
        if (mask.any()):
            ax.scatter(scores[mask, 0], scores[mask, 1], c = colors_np[mask].tolist(), s = 35,
                       marker = shape, edgecolors = "black", linewidths = 0.5)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    ax.set_title(f"{region}: ({nsites} variable sites)")
    ax.legend(handles = handles, frameon = False, loc = "center left", bbox_to_anchor = (1.02, 0.5))

    if (annPCA):
        for i, sample in enumerate(samples_ls):
            ax.annotate(sample, (scores[i, 0], scores[i, 1]), xytext = (4, 4),
                        textcoords = "offset points", fontsize = 6)

    outfile = region.replace(':', '_') + ".png"
    fig.savefig(outfile, dpi = 300, bbox_inches = "tight")
    plt.close(fig)

    return 0

def run_pca() -> int:
    """center the genotype matrix and use an SVD to get PC1 & PC2"""

    global samples, minVar

    # for consistency
    samples_ls = list(samples.keys())

    # create a matrix
    mat = np.array([samples[sample].variants for sample in samples_ls], dtype=float)

    # remove invariant sites
    mat    = mat[:, mat.std(axis = 0) > 0]
    nsites = mat.shape[1]

    if (nsites < minVar or mat.shape[0] < minVar):
        sys.exit(f"{nsites} variable sites found. Minimum required: {minVar}. Exiting..")

    mat           = mat - mat.mean(axis = 0)
    U, S, Vt      = np.linalg.svd(mat, full_matrices=False)
    scores        = U * S
    var_explained = (S**2) / np.sum(S ** 2)

    # now plot
    plot_pca(samples_ls, scores[:, :2], var_explained[:2], nsites)

    return (scores[:, :2], var_explained[:2], nsites)

def main() -> int:
    """entry point to this initiate the entire process"""

    # get the inputs
    parse_command_line()

    # get the sample information
    load_samples()

    # now collect variant positions
    parse_vcf()

    # now do a PCA on each gene
    run_pca()

    return 0
        
if __name__ == '__main__':
    main()
