#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
import copy
from typing import TextIO
from collections import defaultdict

ksize    = 17
aln_file = ''
pop_file = ''
popmap   = dict()

def get_arguments() -> int:
    """grab the input files"""

    global aln_file, pop_file, ksize
    
    desc   = "Helper script to find k-mers where the populations are different in a fasta alignment for each provided sample"
    parser = argparse.ArgumentParser(description=desc)

    # help messages & description
    fhelp = "Fasta alignment file [requierd]"
    phelp = "Population map in tabular format [required]"
    khelp = "K-mer size [default: 17]"

    parser.add_argument("-f", "--fasta", required=True, type=str, help=fhelp)
    parser.add_argument("-p", "--pops",  required=True, type=str, help=phelp)
    parser.add_argument("-k", "--ksize", default=ksize, type=int, help=khelp)

    args = parser.parse_args()
    aln_file = args.fasta
    pop_file = args.pops
    ksize    = args.ksize

    assert os.path.isfile(aln_file), f"Could not locate {aln_file}"
    assert os.path.isfile(pop_file), f"Could not locate {pop_file}"

    assert ksize > 0, f"--ksize cannot be less than 1"

    return 0

def open_file(infile: str) -> TextIO:
    """helper function to create an input stream"""

    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')

def load_pops() -> int:
    """create a population mapping"""

    global pop_file, popmap

    fh  = open_file(pop_file)
    tmp = set()

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 2, "Population map must be 2 columns only [Sample -> Pop]"
        ind         = fields[0] # individual
        pop         = fields[1].strip()
        popmap[ind] = pop
        tmp.add(pop)

    fh.close()

    if (len(popmap) == 0):
        msg = f"No population information was found in {pop_file}"
        sys.exit(msg)
    elif (len(tmp) == 1):
        msg = "Population map only has 1 population."
        sys.exit(msg)

    print("Loaded", len(popmap), "samples")

    return 0

def load_aln() -> dict[str, str]:
    """load the alignment in single-lined fasta format"""

    global aln_file

    # parse and process
    fh        = open_file(aln_file)
    is_header = False # expectation 1 line for header, 1 line for fasta, not multi-line
    samples   = dict() # indvidual -> alignment str
    sample    = ''

    # load the population info
    load_pops()
    global popmap

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        line      = line.strip()
        is_header = not is_header

        if (is_header and line[0] != '>'):
            msg = "Error: Alignment file should not be a multi-aligned fasta"
            sys.exit(msg)

        if (is_header):
            sample = line[1:]
            if (sample not in popmap):
                msg = f"Error: {sample} not found in population map"
                sys.exit(msg)
        else:
            samples[sample] = line

    fh.close()

    # ensure all samples collected
    for sample in popmap.keys():
        if (sample not in samples):
            msg = f"Error: {sample} not found in alignment file"
            sys.exit(msg)

    # ensure all alignments are the same length
    length = -1

    for sample, aln in samples.items():
        if (length == -1):
            length = len(aln)
        elif (length != len(aln)):
            msg = f"Error: Alignments must be the same length"
            sys.exit(msg)

    return samples

def find_pop_kmers(samples: dict[str, str]) -> dict[str, set[str]]:
    """find kmers that are specific to each population"""

    global popmap, ksize

    sample_kmers = defaultdict(set)

    # step 1: collect each samples k-mers and avoid gaps
    for sample, aln in samples.items():
        aln = aln.replace('-', '')
        for i in range(len(aln) - ksize + 1):
            kmer = aln[i:i+ksize]
            sample_kmers[sample].add(kmer)

    # step 2: find which k-mers are shared across the entire pops
    pop_kmer_cnts = defaultdict(lambda: defaultdict(int))
    pop_cnts      = defaultdict(int)

    for sample, kmers in sample_kmers.items():
        pop = popmap[sample]
        for kmer in kmers:
            pop_kmer_cnts[pop][kmer] += 1
        pop_cnts[pop] += 1

    # step 3: find which kmers can be found for all samples
    pop_kmers = defaultdict(set)
    pops      = list()

    for pop, pop_cnt in pop_cnts.items():
        for kmer, kmer_cnt in pop_kmer_cnts[pop].items():
            if (kmer_cnt == pop_cnt):
                pop_kmers[pop].add(kmer)
        pops.append(pop)

    # step 4: find k-mers not shared across the populations
    pop_specific_kmers = dict()

    for i in range(len(pops)):
        cur_pop    = pops[i]
        cur_kmers  = pop_kmers[cur_pop]
        uniq_kmers = copy.copy(cur_kmers)
        if (len(cur_kmers) == 0):
            continue
        for j in range(len(pops)):
            if (j == i):
                continue
            nxt_pop   = pops[j]
            nxt_kmers = pop_kmers[nxt_pop]
            if (len(nxt_kmers) == 0):
                continue
            shared_kmers = uniq_kmers.intersection(nxt_kmers)
            if (len(shared_kmers) > 0):
                uniq_kmers = uniq_kmers.difference(shared_kmers)
        pop_specific_kmers[cur_pop] = uniq_kmers
        
    # now we can return the kmer map
    return pop_specific_kmers

def write_kmers(pop_kmers: dict[str, set[str]]) -> int:
    """Write a simple list of k-mers for each population"""

    for pop, kmers in pop_kmers:
        fh = open(f"{pop}_kmers.txt", 'w')
        for kmer in kmers:
            fh.write(kmer + '\n')
        fh.close()
        print(f"Wrote {len(kmers)} kmers for {pop}")
    
    return 0

def main() -> int:
    """entry point to this helper script"""    

    # retrieve the arguments
    get_arguments()

    # load the pop info
    load_pops()

    # load the alignment
    samples = load_aln()

    # find k-mers specific to each pop
    pop_kmers = find_pop_kmers(samples)

    # now write out the results
    write_kmers(pop_kmers)

    return 0
        
if __name__ == '__main__':
    main()
