#!/usr/bin/env python3
import argparse
import sys
import os
import gzip
from typing import TextIO

#
# this code is specifically designed to take diploid
# calls in males (specifically, male stickleback)
# and convert them to haploid calls. For sites were
# the allele depth (AD) is equal and where the alleles
# are distinct, we will just call the site missing
# instead of being wrong
#
# In this study, the entire PAR region was masked
# Sumarli, A., Behrens, C., Lucas, G., Olufemi, M. J., Bentzen, P., 
# Bell, A. M., ... & Samuk, K. (2025). Strong but diffuse genetic 
# divergence underlies differentiation in an incipient species of 
# marine stickleback. bioRxiv.
#
# The PAR region coordinates are based on the study that assembled
# the Y chromosome
# Peichel, C. L., McCann, S. R., Ross, J. A., Naftaly, A. F., Urton,
# J. R., Cech, J. N., ... & White, M. A. (2020). Assembly of the 
# threespine stickleback Y chromosome reveals convergent signatures 
# of sex chromosome evolution. Genome biology, 21(1), 177.
#
#

# Based on Peichel et al., (2020)
XCHROM  = "XIX"
PAR_END = 2_500_000

# Globals to be set by the user's input
vcf     = ''
popFile = '' # file of <sample> <sex>
outfile = "out.vcf"
males   = list()
cmd     = list()

class Sample:
    __slots__ = ("idx", "hap_count", "mis_count")
    def __init__(self) -> None:
        self.idx       = 0 # index in vcf header
        self.hap_count = 0 # num of het sites converted to haploid
        self.mis_count = 0 # num of het sites converted to missing
    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""

    global outfile, popFile, vcf, cmd
    
    desc  = "Reduce the male heterozgyous sites to haploid calls in the X chromosome in stickleback"
    ohelp = "Name of output file (default out.vcf)"
    shelp = "A two column tsv file with sample and M/F label to flag male individuals"
    vhelp = "Vcf file containing samples found in --samples"
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-s", "--samples", required=True,  type=str, help=shelp)
    parser.add_argument("-v", "--vcf",     required=True,  type=str, help=vhelp)
    parser.add_argument("-o", "--out",     default=outfile, type=str, help=ohelp)

    args = parser.parse_args()
    assert os.path.isfile(args.samples), f"Could not find {args.samples}"
    assert os.path.isfile(args.vcf),     f"Could not find {args.vcf}"

    popFile = args.samples
    vcf     = args.vcf
    outfile = args.out

    if (outfile.endswith(".gz")):
        print("Gzip output file not support. Please use bgzip afterwards")
        outfile = outfile.replace(".gz", '')

    # save the command to put into the vcf header
    for arg in sys.argv:
        if (os.path.isfile(arg)):
            arg = os.path.realpath(arg)
        cmd.append(arg)

    return 0

def open_file(infile: str, io: str) -> TextIO:
    """helper function to create the IO stream for reading or writing files"""
    # io is either r or w 
    return gzip.open(infile, io + 't') if infile.endswith(".gz") else open(infile, io)

def load_males() -> int:
    """helper function to load the ids for each male"""

    global males, popFile

    fh = open_file(popFile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
            continue
        fields = line.strip().split('\t')
        assert len(fields) == 2, f"Error: Expected two columns in --samples. Offending line:\n{line}"
        sample = fields[0]
        sex    = fields[1].upper()
        assert sex in ('M', 'F'), f"Error: Unexpected sex found in --samples:\n{line}"
        if (sex == 'M'):
            males.append(sample)

    fh.close()

    p = os.path.basename(popFile)
    if (len(males) == 0):
        sys.exit(f"No males found in {p}. Nothing to do..")

    print(f"Loaded {len(males)} male samples from {p}")
    return 0

def parse_vcf() -> int:
    """go through the vcf and haploidize male genotypes in the PAR"""

    global vcf, males, XCHROM, PAR_END, outfile, cmd

    male_idxs = list()
    sites     = 0 # total sites visited in PAR
    het_sites = 0 # total sites in PAR that have >=1 hetero male
    male_cnts = {m : Sample() for m in males}
    in_fh     = open_file(vcf, 'r')
    out_fh    = open_file(outfile, 'w')

    for line in in_fh:
        if (len(line) == 0 or line[0] == '\n'):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields     = line.split('\t')
                fields[-1] = fields[-1].strip() # remove end line char
                for i in range(9, len(fields)):
                    sample = fields[i]
                    if (sample in males):
                        male_cnts[sample].idx = i
                        male_idxs.append((sample, i))
                if (len(male_idxs) == 0):
                    sys.exit(f"ERROR: Did not find male sample IDs in the header of {os.path.basename(vcf)}")
                else:
                    # remove samples not in vcf
                    for m in males:
                        if (male_cnts[m].idx == 0):
                            del male_cnts[m]
                out_fh.write(f"##haploidize_command={' '.join(cmd)}\n")
            out_fh.write(line)
            continue
        head   = line.split('\t', 2)
        chrom  = head[0]
        if (chrom != XCHROM):
            out_fh.write(line)
            continue
        # check if this site preceeds the par 
        pos = int(head[1])
        if (pos < PAR_END):
            out_fh.write(line)
            continue
        sites += 1
        fields = line[:-1].split('\t')
        # expected format: GT:DP:RGQ:GQ:MIN_DP:PL
        fmt    = fields[8].split(':')
        ax_idx = len(fmt) # this will trigger out of index error if not found
        for i in range(len(fmt)):
            if (fmt[i] == "AD"):
                ax_idx = i
                break

        # we will write the output line by replacing the entries in fields
        het_found = False
        for (male, i) in male_idxs:
            sample_info = fields[i]
            subfields   = sample_info.split(':')
            sep         = '/' if '/' in subfields[0] else '|'
            alleles     = subfields[0].split(sep)
            if (len(alleles) != 2): # already haploid
                continue
            a1 = alleles[0]
            a2 = alleles[1]
            if (a1 == '.' and a2 == '.'): # missing
                gt = '.'
            elif (a1 == a2): # homozygous
                gt = a1
            elif (a1 == '.' or a2 == '.'): # one of them is missing
                gt = a1 if (a2 == '.') else a2
            else:
                # see if we can choose soley based on allele depth
                try:
                    adepths = subfields[ax_idx].split(',')
                    depth1  = int(adepths[int(a1)])
                    depth2  = int(adepths[int(a2)])
                except (IndexError, ValueError):
                    depth1  = 0
                    depth2  = 0 
                if (depth1 > depth2): # choose a1
                    gt = a1
                elif (depth1 < depth2): # choose a2
                    gt = a2
                else: # set as missing, we can't be sure
                    gt = '.'
                if (gt == '.'):
                    male_cnts[male].mis_count += 1
                else:
                    male_cnts[male].hap_count += 1
                    het_found = True
            subfields[0] = gt
            fields[i]    = ':'.join(subfields)

        het_sites += 1 if het_found else 0
        out_fh.write('\t'.join(fields) + '\n')
            
    in_fh.close()
    out_fh.close()

    if (sites == 0):
        sys.exit("No sites encounterd in PAR region")

    print(f"Total number of sites encountered past the PAR: {sites}")
    print(f"Sites with at least one male het converted to haploid: {het_sites} ({100 * het_sites / sites:.2f}%)")
    print("sample\thets\t%het_of_sites\t%to_haploid\t%to_missing")
    for (male, _) in male_idxs:
        hapCnt = male_cnts[male].hap_count
        misCnt = male_cnts[male].mis_count
        hets   = hapCnt + misCnt 
        d    = max(hets, 1) # avoid dividing by zero
        print(f"{male}\t{hets}\t{100 * hets / sites:.2f}\t{100 * hapCnt / d:.2f}\t{100 * misCnt / d:.2f}")

    return 0

def main() -> int:
    """entry point to this initiate the entire process"""

    # get the inputs
    parse_command_line()

    # load the male samples
    load_males()

    # now perform the haploid calling
    parse_vcf()
   
    return 0
        
if __name__ == '__main__':
    main()
