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
ann         = ''
popmap_file = ''
pops        = list()
popmap      = dict()
genes       = defaultdict(list)

def get_arguments() -> int:
    """function to get the users input"""

    global vcf, popmap_file, ann, region

    # help messages & description
    desc  = "Plot the allele frequency across a specified region for a diploid species"
    rhelp = "region to plot. Format expected CHR:Start-END"
    vhelp = "vcf file"
    ahelp = "Annotation file in gff3 or gtf format"
    phelp = "population map for samples in vcf"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-r", "--region", required=True, type=str, help=rhelp)
    parser.add_argument("-v", "--vcf",    required=True, type=str, help=vhelp)
    parser.add_argument("-p", "--popmap", required=True, type=str, help=phelp)
    parser.add_argument("-a", "--ann",    type=str, help=ahelp, default=ann)

    args = parser.parse_args()

    assert os.path.isfile(args.vcf),    f"Could not locate {args.vcf}"
    assert os.path.isfile(args.popmap), f"Could not locate {args.popmap}"
    assert ':' in args.region,          f"Region expected format: CHR:Start-END"

    if (args.ann != ''):
        assert os.path.isfile(args.ann), f"Could not locate {args.ann}"

    region      = args.region
    vcf         = args.vcf
    popmap_file = args.popmap
    ann         = args.ann

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

def make_attribute_map(attrb: str) -> dict[str, str]:

    amap   = dict()
    attrb  = attrb.strip(' "\n') # remove new line char
    fields = attrb.split(';')

    for field in fields:
        field = field.strip(' "')
        if (field == ''):
            continue
        idx = field.find(' ') # find index of first space
        if (idx == -1 or idx == len(field) - 1):
            idx = field.find('=')
            if (idx == -1 or idx == len(field) - 1):
                continue
        key = field[:idx].lower()
        key = key.strip(' "')
        val = field[idx + 1:]
        val = val.strip(' "')
        amap[key] = val

    return amap

def parse_ann() -> int:
    """parse the annotation file if present"""

    global ann, genes

    if (ann == ''):
        return 1 # nothing to do here

    fh         = gzip.open(ann, "rt") if ann.endswith(".gz") else open(ann, 'r')
    target_reg = parse_region()
    target_chr = target_reg[0] + '\t'
    left       = target_reg[1]
    right      = target_reg[2]

    for line in fh:
        if (len(line) == 0 or line[0] == '#' or line.startswith(target_chr) == False):
            continue
        fields = line.split('\t')
        if (fields[2] != "exon"):
            continue
        start = int(fields[3])
        end   = int(fields[4]) # check for some overlap
        if (start <= right and end >= left):
            attributes = make_attribute_map(fields[8].strip())
            gene_name  = attributes["gene_id"]
            start      = max(start, left)
            end        = min(end, right)
            genes[gene_name].append((start, end))

    fh.close()

    print(f"Found {len(genes)} genes in the target region")

    return 0

def plot_region(freq_map : dict[str, Population]) -> int:
    """Plot each population's allele frequency across the genome"""

    global pops, genes

    # create a shared drawing surface
    target_reg = parse_region()
    npops      = len(pops)
    add_genes  = len(genes) > 0
    nrows      = npops + 1 if add_genes else npops
    ratios     = [1.0] * npops + ([0.7] if add_genes else [])
    fig, axes  = plt.subplots(nrows = nrows, ncols = 1, figsize = (20, 1.5 * npops + 1), 
                            sharex = True, squeeze = False, layout = "constrained",
                            gridspec_kw = {"height_ratios": ratios})

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
            axe.set_title(f"{target_reg[0]}:{target_reg[1]:,}-{target_reg[2]:,}")
            axe.legend(loc = "upper right", markerscale = 3)

        # last row will have the chrom length
        if (add_genes == False and i == npops - 1):
            axe.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / div:,.2f}".rstrip('0').rstrip('.')))
            axe.set_xlabel(f"Position on {target_reg[0]} ({unit})")

    if (add_genes):
        axe    = axes[-1][0]
        colors = ["#3B3B3B", "#2A9D8F", "#E9A23B"]  
        count  = 0
        for gene, exons in genes.items():
            color  = colors[count % len(colors)]
            count += 1
            spans  = [(start, end - start + 1) for start, end in exons]
            axe.broken_barh(spans, (0, 1), facecolors = color, edgecolors = color, linewidth = 0.5)

            # now add the gene label
            # label the first (leftmost) exon, staggered between two heights
            first_exon = min(start for start, _ in exons)
            y = -0.15 if count % 2 == 0 else -0.85
            axe.text(first_exon, y, gene, ha = "left", va = "top",
                     fontsize = 8, color = color, clip_on = True)

        axe.set_ylim(-1.6, 1)
        axe.set_yticks([])
        axe.set_ylabel("Genes", rotation = 90, va = "center")
        axe.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / div:,.2f}".rstrip('0').rstrip('.')))
        axe.set_xlabel(f"Position on {target_reg[0]} ({unit})")


    fig.supylabel("Allele frequency", x = 0.01) # add a shared y-label
    fig.get_layout_engine().set(rect = (0.03, 0, 0.97, 1)) 
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
    else:
        print(f"Found {sites} variant positions in {target_chr}:{left}-{right}")

    return freq_map

def main() -> int:
    """entry point to this application"""

    # get the input files
    get_arguments()

    # grab pop info
    load_population()

    # parse the vcf
    freq_map = parse_vcf()

    # see if any genes will be added to the fig
    parse_ann()

    # now plot
    plot_region(freq_map)
    
    return 0

if __name__ == "__main__":
    main()
