#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
from typing import TextIO

targets  = list()
aln_file = ''
pop_file = ''
popmap   = dict()

def parse_targets(samples: str) -> int:
    """parse the input string to get the samples ids"""

    global targets

    split = samples.split(',')

    for sample in split:
        targets.append(sample.strip(" \'\""))

    if (len(targets) == 0):
        msg = "No targets found"
        sys.exit(msg)
    else:
        print("Targets load:", ", ".join(targets))

    return 0

def get_arguments() -> int:
    """grab the input files"""

    global aln_file, pop_file, targets
    
    desc   = "Helper script to find the positions where the populations are different in a fasta alignment for each provided sample"
    parser = argparse.ArgumentParser(description=desc)

    # help messages & description
    fhelp = "Fasta alignment file [requierd]"
    phelp = "Population map in tabular format [required]"
    thelp = "Comma separated list of samples to find location of fixed differences (e.g., sample1,sample2,etc...) [required]"

    parser.add_argument("-f", "--fasta",   required=True, type=str, help=fhelp)
    parser.add_argument("-p", "--pops",    required=True, type=str, help=phelp)
    parser.add_argument("-t", "--targets", required=True, type=str, help=thelp)

    args = parser.parse_args()
    aln_file = args.fasta
    pop_file = args.pops
    samples  = args.targets

    assert os.path.isfile(aln_file), f"Could not locate {aln_file}"
    assert os.path.isfile(pop_file), f"Could not locate {pop_file}"

    parse_targets(samples)

    return 0

def open_file(infile: str) -> TextIO:
    """helper function to create an input stream"""

    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')

def load_pops() -> int:
    """create a population mapping"""

    global pop_file, popmap, targets

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
        for target in targets:
            if (target not in popmap):
                msg = f"Error: Did not find {target} in population map"
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

def find_fixed_sites(samples: dict[str, str]) -> set[int]:
    """find the fixed differences between the established population"""

    global popmap

    pop_algns = dict() # str -> list[str]
    pops      = list()
    length    = -1
    sites     = set()

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
            if (bool(nucs & pop_nucs)):
                fixed = False
            else:
                nucs.update(pop_nucs)
        if (fixed == False or len(nucs) != len(pops)):
            continue
        sites.add(i)
        tot += 1

    print(f"Out of {length} sites, {tot} were found to be fixed")

    return sites


def write_fixed(sites: set[str], samples: dict[str, str]) -> int:
    """helper function to go through each targets alignment & fixed sites relative to its actual sequence"""

    global targets

    # create a mapping of str -> list of tuples (pos, nuc)
    tmap = dict()
    for target in targets:
        tsites  = list()
        align   = samples[target]
        length  = len(align)
        pos     = 1

        for i in range(length):
            if (i in sites and align[i] != '-'):
                tsites.append((pos, align[i]))
            if (align[i] != '-'):
                pos += 1 # this consumes the read
        tmap[target] = tsites

    # now to write the results

    for target, tsites in tmap.items():
        fh = open(f"{target}_fixed_sites.txt", 'w')
        for tsite in tsites:
            pos = tsite[0]
            nuc = tsite[1]
            fh.write(f"{pos}\t{nuc}\n")
        fh.close()
    
    return 0

def main() -> int:
    """entry point to this helper script"""    

    # retrieve the arguments
    get_arguments()

    # load the alignment
    samples = load_aln()

    # now traverse through the alignments
    sites = find_fixed_sites(samples)

    # now write out the results
    write_fixed(sites, samples)

    return 0
        
if __name__ == '__main__':
    main()
