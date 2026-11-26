#!/bin/env python3

import argparse
import os
import gzip
from typing import TextIO

#
# some simple python code to make a tsv file quickly looking up the gene_id
# and if available gene name, for each protein record
#

ann = ''

def get_arguments() -> int:
    """get the arguments"""

    global ann

    d = "Create a tsv file that is easy to inspect to obtain the gene_id of a protein sequence"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--gtf", help="gtf file", required=True)
    args = parser.parse_args()
    ann  = args.gtf

    assert os.path.isfile(ann), f"Could not locate file: {ann}"

    return 0

def open_file(infile: str) -> TextIO:
    """helper function that only contains a conditional"""
    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')

def get_gene_id(attrb: str) -> str:
    """return the gene_id for this record"""

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
            break
            
    if (gene_id == ''):
        gene_id = id_ # assume an ID= was found

    return gene_id

def make_attrb_map(attrb: str) -> dict[str, str]:
    """return the attributes as key value pairs"""

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
        key = field[:idx]
        key = key.strip(' "')
        val = field[idx + 1:]
        val = val.strip(' "')
        amap[key] = val

    return amap

def write_map(gtf: str, gpMap: dict[str: set[str]], gnames: dict[str, str]) -> int:
    """create an output file"""

    if (".gtf" in gtf):
        r = ".gtf.gz" if gtf.endswith(".gz") else ".gtf"
    elif (".gff" in gtf):
        r = ".gff.gz" if gtf.endswith(".gz") else ".gff"
    b = os.path.basename(gtf).replace(r, ".gpMap.tsv")

    ofh = open(b, 'w')

    for gene_id, prot_ids in gpMap.items():
        gene_name = gnames.get(gene_id, '')
        for prot_id in prot_ids:
            ofh.write(f"{gene_id}\t{prot_id}\t{gene_name}\n")

    ofh.close()

    return 0

def make_map() -> int:
    """One single function that executes all the helper functions for parsing & writing"""

    global ann

    fh     = open_file(ann)
    gpMap  = {}
    gnames = {} # gene names

    for line in fh:
        if (len(line) == 0) or (line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[2].lower() == "gene"):
            attribs  = make_attrb_map(fields[8])
            gene_id  = ''
            if ("gene_id" not in attribs):
                if ("ID" not in attribs):
                    msg = f"Failed to find a gene_id/ID for the following line: {line}"
                    print(msg, end='')
                    continue
                else:
                    gene_id = attribs["ID"]
            else:
                gene_id = attribs["gene_id"]
            if ("gene_name" in attribs):
                gnames[gene_id] = attribs["gene_name"]
            continue
        if (fields[2].lower() != "cds"):
            continue
        attribs   = make_attrb_map(fields[8])
        gene_id   = ''
        prot_id   = ''
        tran_id   = ''
        if ("gene_id" not in attribs):
            if ("ID" not in attribs):
                msg = f"Failed to find a gene_id/ID for the following line: {line}"
                print(msg, end='')
                continue
            else:
                gene_id = attribs["ID"]
        else:
            gene_id = attribs["gene_id"]
        
        if ("protein_id" in attribs):
            prot_id = attribs["protein_id"]
            if ("protein_version" in attribs):
                prot_vs = attribs["protein_version"]
                prot_id = prot_id + '.' + prot_vs
        elif ("transcript_id" in attribs):
            tran_id = attribs["transcript_id"]
            if ("transcript_version" in attribs):
                tran_vs = attribs["transcript_version"]
                tran_id = tran_id + '.' + tran_vs
        elif ("Parent" in attribs):
            tran_id = attribs["Parent"]

        if (tran_id == '' and prot_id == ''):
            msg = f"Failed to find a transcript ID at the following line: {line}"
            print(msg, end='')
        if (gene_id not in gpMap):
            gpMap[gene_id] = set()
        if (prot_id == '' and tran_id != ''):
            prot_id = tran_id
        gpMap[gene_id].add(prot_id)

    fh.close()

    write_map(ann, gpMap, gnames)

    return 0

def main() -> int:
    """get the gtf file and create a tsv file to serve as a map"""

    # get the single argument
    get_arguments()

    # parse & write
    make_map()

    return 0

if __name__ == "__main__":
    main()
