#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict

ann     = '' # single annotation file
vcf     = ''
tarFile = '' # target file containg gene_ids
popFile = '' # file of <sample> <hex color>
outdir  = '.'
minVar  = 2
targets = set()
genes   = defaultdict(list) # chrom -> [Gene class]
samples = dict()            # sample  -> Sample

class Sample:

    __slots__ = ("id", "color", "population", "sex")

    def __init__(self, id: str, color: str, population: str, sex: str) -> None:
        self.id         = id
        self.color      = color
        self.population = population
        self.sex        = sex
class Gene:

    __slots__ = ("chrom", "start", "end", "exons", "genotypes", "id")

    def __init__(self, chrom: str, start: int, end: int) -> None:
        self.chrom     = chrom
        self.start     = start
        self.end       = end
        self.id        = ''
        self.exons     = list()
        self.genotypes = defaultdict(list) # sample: list[genotypes]

    def add_exon(self, start: int, end: int) -> int:
        # just add the exon for now, we will resolve later
        self.exons.append((start, end))
        return 0

    def resolve_exons(self) -> int:

        # merge and collapse any overlapping exons

        if (len(self.exons) <= 1):
            return 1 # nothing to resolve

        self.exons.sort(key = lambda e: e[0])
        resolved = [self.exons[0]]

        i = 1
        while (i < len(self.exons)):
            prev_exon = resolved[-1]
            next_exon = self.exons[i]
            if (next_exon[0] <= prev_exon[1]):
                if (next_exon[1] > prev_exon[1]): # merge
                    resolved[-1] = (prev_exon[0], next_exon[1]) # tuples are immutable
            else:
                resolved.append(next_exon)
            i += 1
        
        self.exons = resolved

        return 0

    def is_exonic(self, pos: int) -> bool:

        for (start, end) in self.exons:
            if (start <= pos <= end):
                return True

        return False

    def add_genotype(self, sample: str, genotype: int) -> int:
        self.genotypes[sample].append(genotype)
        return 0
    
    def has_variants(self) -> bool:
        return len(self.genotypes) > 0
    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""

    global outdir, ann, popFile, vcf, tarFile, minVar, legFile
    
    desc  = "Generate a gene-specific PCA for the target genes using only exonic variants"
    ahelp = "Gene annotation file in GFF3/GTF format"
    ohelp = "Path to output directory to place images"
    ghelp = "Single column list of gene_ids found in annotation"
    shelp = "Four column tsv file containing sample and hex value for PCA color point, population ID, and sex [M/F]"
    lhelp = "Two column tsv file containing hex value and pop ID for PCA Legend"
    vhelp = "Vcf file containing samples found in --samples"
    mhelp = "Minimum number of variant sites required per gene"
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-s", "--samples", required=True,  type=str, help=shelp)
    parser.add_argument("-l", "--legend",  required=True,  type=str, help=lhelp)
    parser.add_argument("-v", "--vcf",     required=True,  type=str, help=vhelp)
    parser.add_argument("-g", "--genes",   required=True,  type=str, help=ghelp)
    parser.add_argument("-a", "--ann",     required=True,  type=str, help=ahelp)
    parser.add_argument("-o", "--outdir",  default=outdir, type=str, help=ohelp)
    parser.add_argument("-m", "--min",     default=minVar, type=int, help=mhelp)

    args = parser.parse_args()
    assert os.path.isfile(args.samples), f"Could not find {args.samples}"
    assert os.path.isfile(args.legend),  f"Could not find {args.legend}"
    assert os.path.isfile(args.vcf),     f"Could not find {args.vcf}"
    assert os.path.isfile(args.genes),   f"Could not find {args.genes}"
    assert os.path.isfile(args.ann),     f"Could not find {args.ann}"
    assert os.path.isdir(args.outdir),   f"Could not find {args.outdir}"
    assert args.min > 1,                 f"--min must be at least 2"

    popFile = args.samples
    legFile = args.legend
    vcf     = args.vcf
    tarFile = args.genes
    ann     = args.ann
    outdir  = args.outdir.rstrip('/')
    minVar  = args.min

    return 0

def load_target_genes() -> int:
    """helper function to load the target gene_ids"""

    global tarFile, targets

    fh = gzip.open(tarFile, "rt") if tarFile.endswith(".gz") else open(tarFile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        gene_id = line.strip()
        targets.add(gene_id)

    fh.close()

    if (len(targets) == 0):
        sys.exit(f"No target gene_ids found in {tarFile}")

    print(f"Loaded {len(targets)} target genes from {tarFile}")
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

def get_gene_id(attrb: str, is_gff: bool) -> str:

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
            
    if (gene_id == '' or is_gff):
        gene_id = id_ # assume an ID= was found

    return gene_id

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
    """pull out the exonic regions for the target genes"""

    global genes, ann, targets

    # now to loop through and parse the annotation
    fh       = gzip.open(ann, "rt") if ann.endswith(".gz") else open(ann, 'r')
    is_gff   = ".gff" in ann
    gene_map = dict()
    total    = 0
    lineNum  = 0
    found    = False
    gene_id  = ''

    for line in fh:
        lineNum += 1
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')

        if (len(fields) != 9):
            cnt = len(fields)
            msg = f"Expected 9 columns per annotation entry. Found {cnt} columns in " + \
                  f"line {lineNum} in {ann}"
            sys.exit(msg) 
        
        feat = fields[2].lower()

        if (feat == "gene"):
            gene_id = get_gene_id(fields[8], is_gff)
            if (gene_id == ''):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann}"
                sys.exit(msg)
            if (gene_id not in targets):
                found = False
                continue
            found          = True
            gene_map[gene_id] = Gene(fields[0], int(fields[3]), int(fields[4]))
            total         += 1
            continue
        elif (feat != "exon" or found == False):
            continue
        
        attrbMap = make_attribute_map(fields[8])

        gene_id = ''
        key1    = "gene_id"
        key2    = "id"
        if (is_gff):
            # the prioritizaiton order should be switched
            key1, key2 = key2, key1
        if (key1 in attrbMap):
            gene_id = attrbMap[key1]
        elif (key2 in attrbMap):
            gene_id = attrbMap[key2]
        if (gene_id == ''):
            msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann}"
            sys.exit(msg)
        start = int(fields[3])
        end   = int(fields[4])
        gene_map[gene_id].add_exon(start, end)         

            
    fh.close()

    if (len(gene_map) == 0):
        sys.exit(f"Did not find any target genes in {ann}")

    print(f"Found {len(gene_map)} out of {len(targets)} genes in {ann}")

    # place each gene into their chromosome buckets
    for gene_id, gene in gene_map.items():
        gene.id = gene_id    # assign ID now
        gene.resolve_exons() # collapse exons
        genes[gene.chrom].append(gene)

    # now sort for binary search later on
    for genes_list in genes.values():
        genes_list.sort(key = lambda g: g.start)
  
    return 0

def get_overlapping_genes(left: int, pos: int, gene_list: list[Gene], gene_ends: list[int], entries: list[Gene]) -> int:
    """for every overlapping gene, add the it to entries list if pos is exonic"""

    while (left < len(gene_list)):
        gene = gene_list[left]
        if (gene.start <= pos <= gene.end):
            if (gene.is_exonic(pos)):
                entries.append(gene)
        elif (gene.start > pos):
            break
        left += 1

    return 0

def parse_vcf() -> int:
    """go through the vcf and find exonic variants for the target genes"""

    global vcf, genes, samples

    samples_ls  = list()
    gene_ends   = list()
    chrom_genes = list()
    curChrom    = ''
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
        head   = line.split('\t', 2)
        chrom  = head[0]
        if (chrom not in genes):
            continue # not target genes on this chrom
        if (curChrom != chrom):
            # we need to update the gene ends
            curChrom    = chrom
            chrom_genes = genes[chrom]
            gene_ends   = [0] * len(chrom_genes)
            longest     = 0
            for i in range(len(chrom_genes)):
                longest      = max(longest, chrom_genes[i].end)
                gene_ends[i] = longest
        # check if this site overlaps any genes
        pos   = int(head[1])
        left  = 0
        right = len(chrom_genes)

        while (left < right):
            mid = left + (right - left) // 2
            if (gene_ends[mid] >= pos):
                right = mid
            else:
                left = mid + 1

        # check if we stopped anywhere
        found = False
        if (left < len(chrom_genes)):
            gene  = chrom_genes[left]
            found = gene.start <= pos
        if (found):
            entries = list() # collect genes where the pos is exonic
            get_overlapping_genes(left, pos, chrom_genes, gene_ends, entries)

            # position was intronic
            if (len(entries) == 0):
                continue
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
                    for j in range(len(entries)):
                        entries[j].add_genotype(sample, score)
                else:
                    missing.append(sample)
            # treat missing value as averages
            score = 0 if total == 0 else total / count
            for sample in missing:
                for j in range(len(entries)):
                    entries[j].add_genotype(sample, score)

    fh.close()

    if (sites == 0):
        sys.exit("No variant sites found in any target genes")
    else:
        print(f"Total number of sites found: {sites}")

    return 0

def run_pca(gene: Gene, samples_ls: list[str]) -> tuple[np.ndarray, np.ndarray, int]:
    """center the genotype matrix and use an SVD to get PC1 & PC2"""

    global minVar

    # create a matrix
    mat = np.array([gene.genotypes[sample] for sample in samples_ls], dtype=float)

    # remove invariant sites
    mat    = mat[:, mat.std(axis = 0) > 0]
    nsites = mat.shape[1]

    if (nsites < minVar or mat.shape[0] < 3):
        return tuple([])

    mat           = mat - mat.mean(axis = 0)
    U, S, Vt      = np.linalg.svd(mat, full_matrices=False)
    scores        = U * S
    var_explained = (S**2) / np.sum(S ** 2)

    return (scores[:, :2], var_explained[:2], nsites)

def plot_pca(gene: Gene, samples_ls: list[str], scores: np.ndarray, var: np.ndarray, nsites: int) -> int:
    """scatter plot of PC1 and PC2 colored by the user's sample colors"""

    global samples, outdir

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

    name    = gene.id.replace("/", '_')              # avoid adding to path
    outpath = f"{outdir}/{name}.pca.png"
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
    ax.set_title(f"{gene.id} ({nsites} variable exonic sites)")
    ax.legend(handles = handles, frameon = False, loc = "center left", bbox_to_anchor = (1.02, 0.5))

    fig.savefig(outpath, dpi = 300, bbox_inches = "tight")
    plt.close(fig)

    return 0

def run_pcas() -> int:
    """run and plot a PCA for each target gene that has enough variants"""

    global genes, samples

    samples_ls = list(samples.keys())
    done       = 0
    skipped    = 0

    for chrom_genes in genes.values():
        for gene in chrom_genes:
            if (gene.has_variants() == False):
                skipped += 1
                continue
            results = run_pca(gene, samples_ls)
            if (len(results) == 0):
                skipped += 1
                continue
            scores = results[0]
            var    = results[1]
            nsites = results[2]
            plot_pca(gene, samples_ls, scores, var, nsites)
            done  += 1

    if (done == 0):
        sys.exit("No genes has enough variable sites for a PCA")

    print(f"Plotted PCA for {done} genes. Skipped {skipped} genes with too few or no variable sites")

    return 0

def main() -> int:
    """entry point to this initiate the entire process"""

    # get the inputs
    parse_command_line()

    # load the target genes
    load_target_genes()

    # get the sample information
    load_samples()

    # construct the gene objects
    parse_ann()

    # now collect variant positions
    parse_vcf()

    # now do a PCA on each gene
    run_pcas()

    return 0
        
if __name__ == '__main__':
    main()
