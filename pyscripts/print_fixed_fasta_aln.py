#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
from typing import TextIO

aln_file = ''
pop_file = ''
popmap   = dict()

def get_arguments() -> int:
    """grab the input files"""

    global aln_file, pop_file
    
    desc   = "Helper script to plot the positions where the populations are different in a fasta alignment"
    parser = argparse.ArgumentParser(description=desc)

    # help messages & description
    fhelp = "Fasta alignment file [requierd]"
    phelp = "Population map in tabular format [required]"

    parser.add_argument("-f", "--fasta", required=True, type=str, help=fhelp)
    parser.add_argument("-p", "--pops",  required=True, type=str, help=phelp)

    args = parser.parse_args()
    aln_file = args.fasta
    pop_file = args.pops

    assert os.path.isfile(aln_file), f"Could not locate {aln_file}"
    assert os.path.isfile(pop_file), f"Could not locate {pop_file}"

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

    else:
        print("Loaded", len(popmap), "samples")

    return 0

def load_aln() -> dict[str, str]:
    """extract the exonic sequences"""

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

def print_fixed_sites(samples: dict[str, str]) -> int:
    """bucket samples to populations and then print"""

    global popmap

    pop_algns = dict() # str -> list[str]
    pops      = list()
    length    = -1

    for sample, aln in samples.items():
        pop = popmap[sample]
        if (pop not in pop_algns):
            pop_algns[pop] = list()
            pops.append(pop)
        pop_algns[pop].append(aln)

        if (length == -1):
            length = len(aln)

    # now iterate over each population
    tot = 0
    for i in range(length):
        fixed = True
        nucs  = set()
        for pop in pops:
            pop_nucs = set()
            alns     = pop_algns[pop]
            for aln in alns:
                pop_nucs.add(aln[i])
                if (len(pop_nucs) > 1):
                    fixed = False
                    break
            if (fixed == False):
                break
            nucs.update(pop_nucs)
        if (fixed == False or len(nucs) != len(pops)):
            continue
        nucs_str = ','.join([n for n in nucs])
        msg      = f"Fixed site {i}: {nucs_str}"
        tot     += 1
        print(msg)

    print(f"Out of {length} sites, {tot} were found to be fixed")

    return 0

def main() -> int:
    """entry point to this helper script"""    

    # retrieve the arguments
    get_arguments()

    # load the alignment
    samples = load_aln()

    # now traverse through the alignments
    print_fixed_sites(samples)

    return 0
        
if __name__ == '__main__':
    main()
