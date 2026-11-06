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
fastas   = list()

def get_arguments() -> int:
    """grab the input files"""

    global aln_file, pop_file, ksize, fastas
    
    desc   = "Helper script to scan for population-specific k-mers obtained from a fasta alignment across fasta sequence files"
    parser = argparse.ArgumentParser(description=desc)

    # help messages & description
    ahelp = "Fasta alignment file [requierd]"
    phelp = "Population map in tabular format [required]"
    khelp = "K-mer size [default: 17]"
    fhelp = "Fasta files (e.g., fa1 fa2 etc..) [required]"

    parser.add_argument("-a", "--algn",   required=True, type=str, help=ahelp)
    parser.add_argument("-p", "--pops",   required=True, type=str, help=phelp)
    parser.add_argument("-f", "--fastas", required=True, type=str, help=fhelp, nargs='*')
    parser.add_argument("-k", "--ksize",  default=ksize, type=int, help=khelp)

    args = parser.parse_args()
    aln_file = args.algn
    pop_file = args.pops
    ksize    = args.ksize
    temp     = args.fastas

    assert os.path.isfile(aln_file), f"Could not locate {aln_file}"
    assert os.path.isfile(pop_file), f"Could not locate {pop_file}"

    assert len(temp) > 0, "--fastas requires input"
    assert     ksize > 0, f"--ksize cannot be less than 1"

    for fa in temp:
        assert os.path.isfile(fa), f"Could not locate {fa}"
        fastas.append(fa)

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
    pop_all_kmers = defaultdict(set) # every k-mer found across each pop, even non-fixed ones

    for sample, kmers in sample_kmers.items():
        pop = popmap[sample]
        for kmer in kmers:
            pop_kmer_cnts[pop][kmer] += 1
        pop_all_kmers[pop].update(kmers)
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

    for cur_pop in pops:
        other_kmers = set()
        for nxt_pop in pops:
            if (cur_pop == nxt_pop):
                continue
            other_kmers.update(pop_all_kmers[nxt_pop])
        pop_specific_kmers[cur_pop] = pop_kmers[cur_pop].difference(other_kmers)
        
    for p, kmers in pop_specific_kmers.items():
        print(f"Found {len(kmers)} kmers unique to {p}")
        
    # now we can return the kmer map
    return pop_specific_kmers

def process_seq(seq: str, pop_kmers: dict[str, set[str]]) -> dict[str, int]:
    """returns counts of either 1 or 0 for each population"""

    global ksize

    counts = {}
    found  = set()

    # generate k-mers for this sequence
    for i in range(len(seq) - ksize + 1):
        kmer = seq[i:i+ksize]
        for pop, kmers in pop_kmers.items():
            if (kmer in kmers):
                found.add(pop)
        # shared by both
        if (len(found) == len(pop_kmers)):
            break

    # try reverse compliment
    rev_trans = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
    seq       = seq.translate(rev_trans)[::-1]
    for i in range(len(seq) - ksize + 1):
            kmer = seq[i:i+ksize]
            for pop, kmers in pop_kmers.items():
                if (kmer in kmers):
                    found.add(pop)
            # shared by both
            if (len(found) == len(pop_kmers)):
                break

    for pop in found:
        counts[pop] = 1


    return counts

def scan_fastas(pop_kmers: dict[str, set[str]]) -> int:
    """print to the console the num of reads each fasta for each category of pop specific k-mers"""

    global fastas

    for fasta in fastas:
        fh   = open_file(fasta)
        cnts = {p : 0 for p in pop_kmers}
        seq  = list()

        # assign a category for both
        cnts["both"] = 0

        for line in fh:
            if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
                continue
            if (line[0] == '>'):
                # new record
                if (len(seq) > 0):
                    seq_cnts = process_seq(''.join(seq), pop_kmers)
                    if (len(seq_cnts) == len(pop_kmers)):
                        cnts["both"] += 1
                    else:
                        for p, cnt in seq_cnts.items():
                            cnts[p] += cnt
                seq.clear()
            else:
                seq.append(line.strip())
        fh.close()

        # process last seq
        if (len(seq) > 0):
            seq_cnts = process_seq(''.join(seq), pop_kmers)
            if (len(seq_cnts) == len(pop_kmers)):
                cnts["both"] += 1
            else:
                for p, cnt in seq_cnts.items():
                    cnts[p] += cnt
            seq.clear()

        print(f"Counts for {os.path.basename(fasta)}:")
        for p, cnt in cnts.items():
            print(f"\tNumber of reads with k-mers shared with {p}:{cnt}")

    return 0

def main() -> int:
    """entry point to this helper script"""    

    # retrieve the arguments
    get_arguments()

    # load the alignment
    samples = load_aln()

    # find k-mers specific to each pop
    pop_kmers = find_pop_kmers(samples)

    # now scan the fasta files for these k-mers
    scan_fastas(pop_kmers)

    return 0
        
if __name__ == '__main__':
    main()
