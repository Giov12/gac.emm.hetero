#!/bin/env python3

import argparse
import os
import gzip
import random
from collections import defaultdict

outdir  = ''
popmap  = ''
pops    = defaultdict(int) # population -> number of samples
samples = list()
orig    = list()           # list of populations following samples' order
perms   = 10_000

def get_arguments() -> int:
    """function to get the users input"""

    global outdir, popmap, perms

    # help messages & description
    desc  = "Generate population maps where samples are randomly assigned to populations"
    phelp = "A two column tsv file containing each sample and its assigned population"
    ohelp = "path to output directory to place population maps"
    nhelp = "Number of permutations to create"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-p", "--popmap", required=True, type=str, help=phelp)
    parser.add_argument("-n", "--number", type=int, help=nhelp)
    parser.add_argument("-o", "--outdir", type=str, help=ohelp, default=outdir)

    args = parser.parse_args()
    assert os.path.isfile(args.popmap), f"Could not find {args.popmap}"
    assert os.path.isdir(args.outdir),  f"Could not find {args.outdir}"
    assert args.number > 0, f"{args.number} must be greater than 0"

    popmap = args.popmap
    outdir = args.outdir.rstrip('/')
    perms  = args.number

    return 0

def load_pop_info() -> int:
    """read in the population may and tally the number of individuals per population"""

    global popmap, pops, samples, orig

    fh = gzip.open(popmap, "rt") if popmap.endswith(".gz") else open(popmap, 'r')

    for line in fh:
        if (len(line) == 0):
            continue
        fields = line.split('\t')
        assert len(fields) == 2, f"Error: Expected two columns in {popmap}\nOffending line: {line}"
        sample = fields[0]
        pop    = fields[1].strip()
        samples.append(sample)
        orig.append(pop)
        pops[pop] += 1

    fh.close()

    print(f"Loaded {len(samples)} samples across {len(pops)}")

    return 0

def permutate() -> int:
    """work horse for this application"""

    global perms, pops, samples, outdir, orig

    # create a set of permutations we already have seen
    seen   = set(tuple(orig))
    labels = orig[:] # create a new copy
    number = len(labels)
    count  = 0

    while (count < perms):
        random.shuffle(labels)
        key = tuple(labels)

        if (key in seen):
            continue # already seen this assignment
        
        seen.add(key)
        out    = f"{outdir}/popmap_{count}.tsv"
        count += 1

        with open(out, 'w') as fh:
            for i in range(number):
                fh.write(f"{samples[i]}\t{labels[i]}\n")

    print(f"Done! Wrote: {perms} permutations in {outdir}/")


    return 0

def main() -> int:
    """entry point to this application"""

    # get the arguments
    get_arguments()

    # load the population information
    load_pop_info()

    # now permutate

    return 0

if __name__ == "__main__":
    main()
