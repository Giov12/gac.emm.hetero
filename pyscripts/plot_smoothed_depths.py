#!/bin/env python3

import argparse
import os
import gzip
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from glob import glob

cov_path    = ''
cov_files   = list()
fai_file    = ''
chroms      = dict() # str -> int
color_map   = dict() # str -> tuple(sample_id, color)
chroms_list = list() # tuples (chrom, int) # same as map
min_len     = 0
xpos_map    = dict() # x-axis position to begin plotting
genomeSize  = 0

def get_arguments() -> int:
    """function to get the users input"""

    global cov_path, fai_file, min_len

    parser = argparse.ArgumentParser(description="Plot the coverages from a set of depth files")
    parser.add_argument("-d", "--dir", required=True, type=str, help="directory containing the depth files")
    parser.add_argument("-f", "--fai", required=True, type=str, help="fasta index file with the lengths of the chromosomes")
    parser.add_argument("-m", "--min", type=int, help="Minimum length of chromosome to include [default: 0]", default=min_len)

    args = parser.parse_args()

    assert os.path.isdir(args.dir),  f"Could not locate {args.dir}"
    assert os.path.isfile(args.fai), f"Could not locate {args.fai}"
    assert args.min >= 0, "--min cannot be less than 0"

    cov_path = args.dir
    fai_file = args.fai
    min_len  = args.min

    return 0

def assign_colors() -> int:
    """assign each file a sample name and color"""

    global cov_files, color_map

    # how many samples do we have?
    count = len(cov_files)

    if (count <= 10):
        cmap   = plt.get_cmap("tab10")
        colors = [cmap(i) for i in range(count)]
    elif (count <= 20):
        cmap   = plt.get_cmap("tab20")
        colors = [cmap(i) for i in range(count)]
    else:
        cmap   = plt.get_cmap("turbo")
        colors = [cmap(i / (count - 1)) for i in range(count)]

    for i, cov_file in enumerate(cov_files):
        fname = os.path.basename(cov_file)
        if (fname.endswith(".gz")):
            sample = fname.replace(".depth.gz", '').replace("Smoothed_", '')
        else:
            sample = fname.replace(".depth", '').replace("Smoothed_", '')
        fields = sample.split('_')
        if (fields[-1][0] == 'S'):
            fields.pop()
        sample = '_'.join(fields)
        color = colors[i]
        color_map[cov_file] = (sample, color)

    return 0

def get_coverages() -> int:
    """get the depth files from the provided path"""

    global cov_path, cov_files

    if (cov_path[-1] == '/'):
        cov_files = glob(f"{cov_path}*.depth*")
    else:
        cov_files = glob(f"{cov_path}/*.depth*")

    assert len(cov_files) > 0, f"Could not locate any .depth files in {cov_path}"

    cov_files.sort() # ensure consistency in plotting

    # collect sample & color info
    assign_colors()

    return 0

def load_chrom_lengths() -> int:
    """read in the chromosome lengths from the index file"""

    global fai_file, chroms, chroms_list, min_len, genomeSize

    fh = open(fai_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        chrom  = fields[0]
        length = int(fields[1])
        if (length >= min_len):
            chroms[chrom] = length
            genomeSize   += length
            chroms_list.append((chrom, length))

    fh.close()

    if (len(chroms) == 0):
        msg = f"No chromosomes found with a minumum length of {min_len}"
        sys.exit(msg)

    # now set the x positions
    chroms_list.sort(key=lambda x: x[1], reverse=True) # largest chroms first
    xstart = 0
    for chrom, chrom_len in chroms_list:
        xpos_map[chrom] = xstart
        xstart         += chrom_len

    return 0


def plot_coverage(infile: str, color: str, ax) -> int:
    """plot a single depth file into the shared ax object"""

    global chroms, xpos_map

    fh        = gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, "r")
    curChrom  = ''
    scores    = list()
    positions = list()
    xpos      = 0
    
    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line[:-1].split('\t') # removes new line character
        chrom  = fields[0]
        if (chrom not in chroms):
            continue # not drawing
        pos   = int(fields[1])
        depth = float(fields[2])
    
        if (curChrom == ''):
            curChrom = chrom
            xpos     = xpos_map[chrom]
        if (curChrom == chrom):
            positions.append(pos + xpos)
            scores.append(depth)
        else:
            ax.plot(positions, scores, color = color, linewidth = 0.5, alpha = 0.8)
            positions.clear()
            scores.clear()
            curChrom = chrom
            xpos     = xpos_map[chrom]
            positions.append(pos + xpos)
            scores.append(depth)
    
    fh.close()
    
    # plot last chrom
    if (len(positions) > 0):
        ax.plot(positions, scores, color = color, linewidth = 0.5, alpha = 0.5)
        positions.clear()
        scores.clear()

    return 0

def plot_coverages() -> int:
    """iteratively plot each depth file"""

    global cov_files, genomeSize, color_map

    # create a shared drawing surface
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (16, 5))

    ax.set_ylim(0, 100)
    ax.set_xlim(0, genomeSize)

    for cov_file in cov_files:
        color = color_map[cov_file][1]
        plot_coverage(cov_file, color, ax)

    # now draw the chromosome & x-axis labels 
    global chroms_list

    current_pos = 0
    xlabels     = list()
    xticks      = list()

    for i, chrom_tuple in enumerate(chroms_list, start=1):
        chrom     = chrom_tuple[0]
        chrom_len = chrom_tuple[1]
        # draw the x-axis chrom label
        xpos = current_pos + (chrom_len / 2)
        xticks.append(xpos)
        xlabels.append(chrom)
        current_pos += chrom_len
        if (i != len(chroms_list)): # not the last chrom
            ax.axvline(current_pos, color = "black", linestyle = "--")

    # add chrom labels
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    # draw the legend
    handles = list()
    labels  = list()
    for sample_info in color_map.values():
        sample = sample_info[0]
        color  = sample_info[1]
        handle = Line2D([0], [0], color = color, linewidth = 2.0)
        labels.append(sample)
        handles.append(handle)
    ncol = min(len(labels), 6)

    ax.legend(handles, labels, loc = "lower center", bbox_to_anchor = (0.5, 1.0),
              ncol = ncol, fancybox = True, shadow = True)

    plt.savefig("genome_coverage.png", dpi = 500, bbox_inches = "tight")

    return 0

def main() -> int:
    """entry point to this application"""

    # get the two inputs
    get_arguments()

    # load the chromosome lengths & get the depth files
    load_chrom_lengths()
    get_coverages()

    # now plot
    plot_coverages()

    return 0

if __name__ == "__main__":
    main()
