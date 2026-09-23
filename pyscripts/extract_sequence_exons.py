#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
import textwrap

target_file = ''
ref_fasta   = ''
ann_file    = ''
outName     = "candidate_seqs.fa.gz"
targets     = set()
chroms      = set()
genes       = dict()

# create once
_DNA_COMPLEMENT = str.maketrans("ATCGatcgNn", "TAGCtagcNn")

class Seq:
    
    __slots__ = ("header", "seq")

    def __init__(self, header: str, seq: str) -> None:
        self.header = header
        self.seq    = seq

    def __len__(self) -> int:
        return len(self.seq)

class Transcript:

    __slots__ = ("id", "exons")

    def __init__(self, id_: str) -> None:
        self.id    = id_
        self.exons = list()

    def add_exons(self, start: int, end: int) -> None:
        self.exons.append((start - 1, end)) # set the indices to be 0-based

    def get_exons(self) -> list[tuple[int, int]]:
        self.exons.sort(key = lambda x: x[0])
        return self.exons
    
class Gene:

    __slots__ = ("id", "chr", "transcripts", "strand")

    def __init__(self, id_: str, chrom: str, strand: str) -> None:
        self.id          = id_
        self.chr         = chrom
        self.transcripts = list()
        self.strand      = strand

    
    def add_transcript(self, transcript_id: str, start: int, end: int) -> int:
        for i in range(len(self.transcripts)):
            if (self.transcripts[i].id == transcript_id):
                self.transcripts[i].add_exons(start, end)
                return 0
            
        #
        # if not executed, add the transcript
        #
        transcript = Transcript(transcript_id)
        transcript.add_exons(start, end)
        self.transcripts.append(transcript)
        return
    
    def get_transcripts(self) -> list[Transcript]:
        return self.transcripts
    

def get_arguments() -> int:
    """grab the input files & optional output file name"""

    global target_file, ref_fasta, outName, ann_file
    
    desc   = "Generate concatenated exonic sequences for genes/transcripts of interest"
    parser = argparse.ArgumentParser(description=desc)

    # help messages & description
    thelp = "A list of gene_ids to extract from the reference genome"
    ahelp = "Annotation file in GFF3/GTF format"
    rhelp = "Reference genome in fasta format to extract exonic regions from"
    ohelp = "Output name [optional]"

    parser.add_argument("-t", "--targets", required=True, type=str, help=thelp)
    parser.add_argument("-r", "--ref",     required=True, type=str, help=rhelp)
    parser.add_argument("-a", "--ann",     required=True, type=str, help=ahelp)
    parser.add_argument("-o", "--out",     type=str,      help=ohelp, default=outName)

    args = parser.parse_args()
    target_file = args.targets
    ref_fasta   = args.ref
    ann_file    = args.ann
    outName     = args.out

    assert os.path.isfile(target_file), f"Could not locate {target_file}"
    assert os.path.isfile(ref_fasta), f"Could not locate {ref_fasta}"
    assert os.path.isfile(ann_file), f"Could not locate {ann_file}"

    return 0

def get_targets() -> int:
    """population the targets container"""

    global target_file, targets

    fh = gzip.open(target_file, "rt") if target_file.endswith(".gz") else open(target_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        line = line.strip()
        targets.add(line)

    fh.close()

    if (len(targets) == 0):
        msg = f"No target genes were provided in {target_file}"
        sys.exit(msg)

    else:
        print("Loaded", len(targets), "target genes")

    return 0

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

    map    = dict()
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
        map[key] = val

    return map

def get_target_exons() -> int:
    """grab the transcripts ids for the target genes"""

    global targets, ann_file, chroms, genes

    # now to loop through and parse the annotation
    fh      = gzip.open(ann_file, "rt") if ann_file.endswith(".gz") else open(ann_file, 'r')
    lineNum = 0
    count   = 0

    for line in fh:
        lineNum += 1
        if (len(line) == 0 or line[0] == '#'):
            continue

        fields = line.split('\t')
        feat   = fields[2].lower()

        if (feat == "gene"):
            gene_id = get_gene_id(fields[8])
            if (gene_id == ''):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann_file}"
                sys.exit(msg)
            if (gene_id in targets):
                genes[gene_id] = Gene(gene_id, fields[0], fields[6])
                chroms.add(fields[0])
            continue

        elif (feat != "exon"):
            continue

        attrbMap = make_attrb_map(fields[8])

        if ("gene_id" not in attrbMap):
            if ("ID" not in attrbMap):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann_file}"
                sys.exit(msg)
            else:
                gene_id = attrbMap["ID"]
        else:
            gene_id = attrbMap["gene_id"]

        if (gene_id not in targets):
            continue
        
        start   = int(fields[3])
        end     = int(fields[4])
        tran_id = ''
        if ("transcript_id" in attrbMap):
            tran_id = attrbMap["transcript_id"]
            if ("transcript_version" in attrbMap):
                tran_vs = attrbMap["transcript_version"]
                tran_id = tran_id + '.' + tran_vs
        if (tran_id == ''):
            msg = f"Failed to find a transcript ID at line {lineNum} in {ann_file}"
            sys.exit(msg)
        genes[gene_id].add_transcript(tran_id, start, end)
        count += 1
            
    fh.close()

    print(f"Loaded {len(genes)} genes from {ann_file}")

    return 0

def get_header_id(header: str) -> str:
    """return the first characters before any spacing"""

    id_ = ''
    idx = header.find(' ')

    if (idx == -1):
        id_ = header[1:]
    else:
        id_ = header[1:idx]

    return id_

def reverse_transcribe(seq: str) -> str:
    """reverse transcribe a sequence"""

    return seq.translate(_DNA_COMPLEMENT)[::-1]

def make_sequences() -> int:
    """extract the exonic sequences"""

    global chroms, ref_fasta, outName, genes

    # parse and process
    fh      = gzip.open(ref_fasta, "rt") if ref_fasta.endswith(".gz") else open(ref_fasta, 'r')
    outfh   = gzip.open(outName, "wt") if outName.endswith(".gz") else open(outName, 'w')
    lineNum = 0
    curRec  = ''
    curSeq  = list()
    count   = 0

    for line in fh:

        lineNum += 1

        if (len(line) == 0 or line[0] == '#'):
            continue
        line = line.strip()

        if (line[0] == '>'):
            if (curRec == ''):
                curRec = get_header_id(line)
                if (curRec not in chroms):
                    curRec = ''
            else:
                seq = ''.join(curSeq)
                for gene_id, gene in genes.items():
                    if (gene.chr != curRec):
                        continue
                    transcripts = gene.get_transcripts()
                    for transcript in transcripts:
                        header = '>' + gene_id + '_' + transcript.id + '\n'
                        exons  = transcript.get_exons()
                        eseq   = list()
                        # exons are tuples of start:end positions
                        for exon in exons:
                            ex = seq[exon[0]:exon[1]]
                            eseq.append(ex)
                        Eseq = ''.join(eseq)
                        if (gene.strand == '-'):
                            Eseq = reverse_transcribe(Eseq)
                        Eseq  = textwrap.fill(Eseq, width=60) + '\n'
                        outfh.write(header)
                        outfh.write(Eseq)
                        count += 1

                # update
                curSeq.clear()
                seq    = ''
                curRec = get_header_id(line)
                if (curRec not in chroms):
                    curRec = ''

        elif (curRec != ''):
            curSeq.append(line)

    fh.close()

    if (curRec != ''):
        seq = ''.join(curSeq)
        for gene_id, gene in genes.items():
            if (gene.chr != curRec):
                continue
            transcripts = gene.get_transcripts()
            for transcript in transcripts:
                header = '>' + gene_id + '_' + transcript.id + '\n'
                exons  = transcript.get_exons()
                eseq   = list()
                # exons are tuples of start:end positions
                for exon in exons:
                    ex = seq[exon[0]:exon[1]]
                    eseq.append(ex)
                Eseq = ''.join(eseq)
                if (gene.strand == '-'):
                    Eseq = reverse_transcribe(Eseq)
                Eseq  = textwrap.fill(Eseq, width=60) + '\n'
                outfh.write(header)
                outfh.write(Eseq)
                count += 1

    outfh.close()

    return 0

def main() -> int:
    """entry point to this helper script"""    

    # retrieve the arguments
    get_arguments()

    # load the target ids
    get_targets()

    # parse the annotation for the genomic coordinates
    get_target_exons()

    # parse the reference genome 
    make_sequences()
    
    return 0
        
if __name__ == '__main__':
    main()
