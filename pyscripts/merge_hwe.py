#!/bin/env python3

import argparse
import os
import gzip

hwe_files = list()
outfile   = "merged_hwe.tsv.gz"

def get_arguments() -> int:
    """function to get the users input"""

    global hwe_files, outfile

    # help messages & description
    desc  = "Merge hwe tables with matching SNPs generated from vcftools --hardy"
    hhelp = "One or more tables generated from vcftools hardy"
    ohelp = "name of output [optional]"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-H", "--hwe-files", required=True, nargs='*', help=hhelp)
    parser.add_argument("-o", "--out", type=str, help=ohelp, default=outfile)

    args = parser.parse_args()
    hwe_files = args.hwe_files
    outfile   = args.out

    assert len(hwe_files) > 0, "--hwe-files requires at least one file"

    for hfile in hwe_files:
        assert os.path.isfile(hfile), f"Could not locate {hfile}"

    if (outfile.endswith(".gz") == False):
        outfile = outfile + ".gz"

    return 0

def merge() -> int:
    """single function that executes the merging"""
    global outfile, hwe_files

    # consistent order
    hwe_files.sort()

    lines = list() # list of lists

    for i, hfile in enumerate(hwe_files):
        fh = gzip.open(hfile, "rt") if hfile.endswith(".gz") else open(hfile, 'r')

        for j, line in enumerate(fh):
            if (i == 0): # first file
                lines.append(list())
            if (j == 0): # column
                if (i == 0):
                    lines[0].append("CHR")
                    lines[0].append("POS")
                # add the sample name
                fname = os.path.basename(hfile).replace(".hwe", '').replace(".gz", '')
                lines[0].append(fname)
            else:
                fields = line.split('\t')
                if (i == 0):
                    lines[j].append(fields[0])
                    lines[j].append(fields[1])
                lines[j].append(fields[2])
            
        fh.close()

    ofh = gzip.open(outfile, "wt")

    for line in lines:
        line = '\t'.join(line) + '\n'
        ofh.write(line)

    ofh.close()
    return 0

def main() -> int:
    """entry point to this application"""

    # get the two inputs
    get_arguments()

    # single function to execute everything
    merge()

    return 0

if __name__ == "__main__":
    main()
