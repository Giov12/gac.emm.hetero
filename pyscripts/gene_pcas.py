#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
from collections import defaultdict

outdir  = '' # output directory
ann     = '' # single annotation file
vcf     = ''
tarFile = '' # target file
popFile = '' # file of <sample> <hex color>
targets = set()
genes   = defaultdict(list) # chrom -> [Gene class]
colmap  = dict()            # sample  -> color
ann     = list()            # annotation file

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
        resolved = list()

        self.exons = resolved

    def is_exonic(self, pos: int) -> bool:

        for (start, end) in self.exons:
            if (start <= pos <= end):
                return True

        return False
    
def assert_file_exists(file_path: str) -> int:
    """avoid copy and pasting the same assert function"""

    assert os.path.isfile(file_path), f"Could not locate {file_path}"

    return 0


def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""

    global outdir, ann, popFile, vcf, tarFile
    
    desc  = "Generate a gene-specific PCA for the target genes using only exonic variants"
    ahelp = "Gene annotation file in GFF3/GTF format"
    ohelp = "Path to output directory to place images"
    ghelp = "Single column list of gene_ids found in annotation"
    shelp = "Two column tsv file containing sample and hex value for PCA color point"
    vhelp = "Vcf file containing samples found in --samples"
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-s", "--samples", required=True, type=str, help=shelp)
    parser.add_argument("-v", "--vcf",     required=True, type=str, help=vhelp)
    parser.add_argument("-g", "--genes",   required=True, type=str, help=ghelp)
    parser.add_argument("-a", "--ann",     required=True, type=str, help=ahelp)
    parser.add_argument("-o", "--outdir",  default=outdir, type=str, help=ohelp)

    args = parser.parse_args()
    assert os.path.isfile(args.samples), f"Could not find {args.samples}"
    assert os.path.isfile(args.vcf),     f"Could not find {args.vcf}"
    assert os.path.isfile(args.genes),   f"Could not find {args.genes}"
    assert os.path.isfile(args.ann),     f"Could not find {args.ann}"
    assert os.path.isdir(args.outdir),   f"Could not find {args.outdir}"

    popFile = args.samples
    vcf     = args.vcf
    tarFile = args.genes
    ann     = args.ann
    outdir  = args.outdir

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

def load_color_map() -> int:
    """load each sample and their color in the pca plots"""

    global colmap, popFile

    fh = gzip.open(popFile, "rt") if popFile.endswith(".gz") else open(popFile, 'r')

    for line in fh:
        if (len(line) == 0):
            continue
        fields = line.split('\t')
        assert len(fields) == 2, f"Malformed line detected in {popFile}:\n{line}"
        assert fields[1][0] == '#', f"Second column in {popFile} should be a hex value. Found {fields[1]}"
        sample         = fields[0]
        color          = fields[1].strip()
        colmap[sample] = color

    fh.close()

    if (len(colmap) == 0):
        sys.exit(f"No samples found in {popFile}")

    print(f"Loaded {len(colmap)} samples from {popFile}")

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
    lineNum  = 0
    found    = False
    gene_id  = ''

    for line in fh:
        lineNum += 1
        if (len(line) == '' or line[0] == '#'):
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

    if (len(genes) == 0):
        sys.exit(f"Did not find any target genes in {ann}")

    print(f"Found {len(genes)} out of {len(targets)} genes in {ann}")

    # place each gene into their chromosome buckets
    for gene_id, gene in gene_map.items():
        gene.id = gene_id # assign ID now
        genes[gene.chrom].append(gene)

    # now sort for binary search later on
    for genes_list in genes.values():
        genes_list.sort(key = lambda g: g.start)
  
    return 0

def get_overlapping_genes(left: int, pos: int, gene_list: list[Gene], gene_ends: list[int], entries: list[Gene]) -> int:
    """for every overlapping gene, add the it to entries list if pos is exonic"""

    return 0

def parse_vcf() -> int:
    """go through the vcf and find exonic variants for the target genes"""

    global vcf, genes, colmap

    nsamples = len(colmap) # expected number of samples

def main() -> int:
    """entry point to this initiate the entire process"""
    
    parse_command_line()


    return 0
        
if __name__ == '__main__':
    main()
