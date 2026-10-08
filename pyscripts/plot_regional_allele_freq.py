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

region      = ''
vcf         = ''
popmap_file = ''
pops        = list()
popmap      = dict()

def get_arguments() -> int:
    """function to get the users input"""

    global vcf, popmap_file, fai_file, region

    # help messages & description
    desc  = "Plot the allele frequency across a specified region for a diploid species"
    rhelp = "region to plot. Format expected CHR:Start-END"
    vhelp = "vcf file"
    phelp = "population map for samples in vcf"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-r", "--region", required=True, type=str, help=rhelp)
    parser.add_argument("-v", "--vcf",    required=True, type=str, help=vhelp)
    parser.add_argument("-p", "--popmap", required=True, type=str, help=phelp)

    args = parser.parse_args()

    assert os.path.isfile(args.vcf),    f"Could not locate {args.vcf}"
    assert os.path.isfile(args.popmap), f"Could not locate {args.popmap}"
    assert ':' in args.region,          f"Region expected format: CHR:Start-END"

    region      = args.region
    vcf         = args.vcf
    popmap_file = args.popmap

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

def plot_region(freq_map : dict[str, Population]) -> int:
    """Plot each population's allele frequency across the genome"""

    global pops

    # create a shared drawing surface
    target_reg = parse_region()
    npops      = len(pops)
    fig, axes  = plt.subplots(nrows = len(pops), ncols = 1, figsize = (20, 1.5 * npops + 1), 
                            sharex = True, squeeze = False, layout = "constrained")

    # estimate the x-axis label scaling
    span = target_reg[2] - target_reg[1]
    if (span >= 5_000_000):
        div, unit = 1e6, "Mb"
    elif (span >= 5_000):
        div, unit = 1e3, "kb"
    else:
        div, unit = 1, "bp"

    for i in range(npops):
        pop = pops[i]
        axe = axes[i][0]
        paf = freq_map[pop] # population allele frequency
        axe.set_ylim(0, 1.0)
        axe.set_xlim(target_reg[1], target_reg[2])
        axe.scatter(paf.positions, paf.allele1_freq, color = "#D6455D", s = 6, label = "REF")
        axe.scatter(paf.positions, paf.allele2_freq, color = "#2F5DA8", s = 6, label = "ALT")

        # add population name on y-axis
        axe.set_ylabel(pop, rotation = 90, va = "center")

        # title will be chromosome
        if (i == 0):
            axe.set_title(target_reg[0])
            axe.legend(loc = "upper right", markerscale = 3)

        # last row will have the chrom length
        if (i == npops - 1):
            axe.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / div:,.2f}".rstrip('0').rstrip('.')))
            axe.set_xlabel(f"Position on {target_reg[0]} ({unit})")

    fig.supylabel("Allele frequency", x = 0.005) # add a shared y-label
    fig.tight_layout()
    fig.savefig(f"{target_reg[0]}_{target_reg[1]}-{target_reg[2]}.png", dpi = 150)
    plt.close(fig)
        
    # free memory
    for pop in freq_map.values():
        pop.clear()

    return 0

def parse_vcf() -> dict[str, Population]:
    """this function will parse the vcf file and return the allele frequencies for each population"""

    global vcf, popmap, pops

    popindices = list()
    freq_map   = {pop : Population() for pop in pops}
    fh         = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')
    target_reg = parse_region()
    target_chr = target_reg[0]
    left       = target_reg[1]
    right      = target_reg[2]
    sites      = 0

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
        
        if (chrom != target_chr):
            if (sites > 0):
                break # already got the sites
            continue

        pos = int(fields[1])

        if ((left <= pos <= right) == False):
            continue

        sites += 1 # count how many records we've encountered

        # calculate population-specific allele frequencies
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
    
    fh.close()

    if (sites == 0):
        sys.exit(f"No records found in {target_chr}:{left}-{right}")

    return freq_map

def main() -> int:
    """entry point to this application"""

    # get the input files
    get_arguments()

    # grab pop info
    load_population()

    # parse the vcf
    freq_map = parse_vcf()

    # now plot
    plot_region(freq_map)
    
    return 0

if __name__ == "__main__":
    main()
