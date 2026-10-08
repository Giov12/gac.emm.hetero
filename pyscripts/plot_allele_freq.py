#!/bin/env python3

import argparse
import os
import gzip
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from collections import defaultdict

class Population:
    def __init__(self) -> None:
        self.positions    = list()
        self.allele1_freq = list()
        self.allele2_freq = list()

    def clear(self) -> None:
        self.positions.clear()
        self.allele1_freq.clear()
        self.allele2_freq.clear()

class Chrom:
    def __init__(self, name: str, length: int) -> None:
        self.name   = name
        self.length = length

vcf         = ''
popmap_file = ''
fai_file    = ''
pops        = list()
chroms      = dict() # str -> Chrom
popmap      = dict()

def get_arguments() -> int:
    """function to get the users input"""

    global vcf, popmap_file, fai_file

    # help messages & description
    desc  = "Plot the allele frequency across the genome for a diploid species"
    vhelp = "vcf file"
    phelp = "population map for samples in vcf"
    fhelp = "fasta index file with the lengths of the chromosomes"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-v", "--vcf",    required=True, type=str, help=vhelp)
    parser.add_argument("-p", "--popmap", required=True, type=str, help=phelp)
    parser.add_argument("-f", "--fai",    required=True, type=str, help=fhelp)

    args = parser.parse_args()

    assert os.path.isfile(args.vcf),    f"Could not locate {args.vcf}"
    assert os.path.isfile(args.popmap), f"Could not locate {args.popmap}"
    assert os.path.isfile(args.fai),    f"Could not locate {args.fai}"

    vcf         = args.vcf
    popmap_file = args.popmap
    fai_file    = args.fai

    return 0

def load_chrom_lengths() -> int:
    """read in the chromosome lengths from the index file"""

    global fai_file, chroms

    fh = open(fai_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields        = line.split('\t')
        chrom         = fields[0]
        length        = int(fields[1])
        chroms[chrom] = Chrom(chrom, length)

    fh.close()

    if (len(chroms) == 0):
        msg = f"No chromosomes found in {fai_file}"
        sys.exit(msg)

    return 0

def load_population() -> int:
    """parse the tsv file containing the population assignments"""

    global popmap_file, popmap, pops

    fh = gzip.open(popmap_file, "rt") if popmap_file.endswith(".gz") else open(popmap_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        parts = line.split('\t')
        assert len(parts) == 2, f"Malformed line found in {popmap_file}: \n{line}"
        sample         = parts[0]
        pop            = parts[1].strip()
        popmap[sample] = pop
        if (len(pops) == 0):
            pops.append(pop)
        elif (pops[-1] != pop):
            pops.append(pop)

    fh.close()

    if (len(popmap) == 0):
        msg = f"No samples found in {popmap_file}"
        sys.exit(msg)
    else:
        print(f"Loaded {len(popmap)} samples across {len(pops)} populations")

    return 0

def plot_chromosome(freq_map : dict[str, Population], chrom: str, pops: list[str]) -> int:
    """Plot each population's allele frequency across the genome"""

    global chroms

    # create a shared drawing surface
    npops     = len(pops)
    fig, axes = plt.subplots(nrows = len(pops), ncols = 1, figsize = (20, 1.5 * npops + 1), 
                            sharex = True, squeeze = False, layout = "constrained")
    chrom_len = chroms[chrom].length

    for i in range(npops):
        pop = pops[i]
        axe = axes[i][0]
        paf = freq_map[pop] # population allele frequency
        axe.set_ylim(0, 1.0)
        axe.set_xlim(0, chrom_len)
        axe.scatter(paf.positions, paf.allele1_freq, color = "#D6455D", s = 3, label = "REF")
        axe.scatter(paf.positions, paf.allele2_freq, color = "#2F5DA8", s = 3, label = "ALT")

        # add population name on y-axis
        axe.set_ylabel(pop, rotation = 90, va = "center")

        # title will be chromosome
        if (i == 0):
            axe.set_title(chrom)
            axe.legend(loc = "upper right", markerscale = 3)

        # last row will have the chrom length
        if (i == npops - 1):
            axe.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:g}"))
            axe.set_xlabel("Position (Mb)")

    fig.supylabel("Allele frequency", x = 0.005) # add a shared y-label
    fig.tight_layout()
    fig.savefig(f"{chrom}.png", dpi = 150)
    plt.close(fig)
        
    # free memory
    for pop in freq_map.values():
        pop.clear()

    return 0

def parse_vcf_and_plot() -> int:
    """this function will parse the vcf file and plot each chromosome"""

    global vcf, popmap, pops

    popindices = list()
    curChrom   = ''
    freq_map   = {pop : Population() for pop in pops}
    fh         = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields = line.strip().split('\t')
                for i in range(9, len(fields)):
                    popindices.append(popmap[fields[i]])
            continue
        fields = line.strip().split('\t')
        chrom  = fields[0]
        if (curChrom == ''):
            curChrom = chrom
        if (curChrom != chrom):
            plot_chromosome(freq_map, curChrom, pops)
            curChrom = chrom
        # calculate population-specific allele frequencies
        pos  = int(fields[1])
        cnts = {p : defaultdict(int) for p in pops}
        for i in range(9, len(fields)):
            pop      = popindices[i - 9]
            genotype = fields[i].split(':')[0]
            allele1  = ''
            allele2  = ''
            k        = 0
            missing  = False
            for j in range(len(genotype)):
                if (genotype[j] == '.'):
                    missing = True
                    break
                if (genotype[j] == '/' or genotype[j] == '|'):
                    k       = j + 1
                    allele1 = genotype[:j]
            if (missing == False):
                allele2 = genotype[k:]
                cnts[pop][allele1] += 1
                cnts[pop][allele2] += 1
        # now add the population frequencies
        for pop, cnt in cnts.items():
            cnt1 = cnt['0']
            cnt2 = cnt['1']
            tot  = cnt1 + cnt2
            if (tot == 0):
                continue # completely missing site
            freq_map[pop].allele1_freq.append(cnt1 / tot)
            freq_map[pop].allele2_freq.append(cnt2 / tot)
            freq_map[pop].positions.append(pos)

    # plot last chrom
    plot_chromosome(freq_map, curChrom, pops)

    fh.close()

    return 0

def main() -> int:
    """entry point to this application"""

    # get the input files
    get_arguments()

    # load the chromosome lengths
    load_chrom_lengths()

    # grab pop info
    load_population()

    # now parse & plot
    parse_vcf_and_plot()
    
    return 0

if __name__ == "__main__":
    main()
