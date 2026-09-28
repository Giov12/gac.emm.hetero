#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <zlib.h>
#include <iomanip>

using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::vector;
using std::cerr;
using std::cout;
using std::find;
using std::stoi;
using std::sort;
using std::unordered_map;
using std::setprecision;


//
// find sites where the population of
// interest is homozygous for a particular
// allele that is not in a homozygous state
// in the other samples
//

typedef unsigned int uint;

enum SnpType {Exonic,Intronic};

struct Exon {
    uint start;
    uint end;
};

struct Attribute {
    string key;
    string val;
};

class Gene {

public:
    string       id;
    string       name;
    uint         start;
    uint         end;
    vector<Exon> exons;

    //
    // empty constructor
    //
    Gene (string id_, string name, uint start, uint end){
        this->id    = id_;
        this->name  = name;
        this->start = start;
        this->end   = end;
    };
    ~Gene(){
        this->exons.clear();
    }

    void add_exon(Exon exon){
        //
        // just add the exon for now
        // & then we will resolve
        // the exons at the end
        //

        this->exons.push_back(exon);
    }

    void resolve_exons(void){
        //
        // de-duplicate & merge overlapping exons
        //

        // handles empty and single exon cases
        if (this->exons.size() <= 1){
            return;
        }

        vector<Exon> resolved;

        const int count = this->exons.size();
        resolved.reserve(count);

        //
        // sort to then just go exon by exon
        //
        sort(this->exons.begin(), this->exons.end(), []
            (const Exon &exon1, const Exon &exon2){
                if (exon1.start == exon2.start){
                    return exon1.end < exon2.end;
                }
                return exon1.start < exon2.start;
            }
        );

        resolved.push_back(exons.front());
        int i = 1;

        while (i < count){
            Exon &prev = resolved.back();
            Exon &next = this->exons[i];

            // is there overlap?
            if (next.start <= prev.end){
                if (next.end > prev.end){ // merge if true
                    prev.end = next.end;
                }
            }
            else {
                resolved.push_back(next);
            }
            i++;
        }
        this->exons = resolved;
    }

    SnpType classify(uint pos){
        //
        // exons are ordered,
        // so just see if there is an overlap
        //

        for (uint i = 0; i < this->exons.size(); i++){
            const Exon &e = this->exons[i];
            if (e.start <= pos && pos <= e.end){
                return SnpType::Exonic;
            }
        }

        // no overlap
        return SnpType::Intronic;
    }
};

struct SNP {
    uint pos;
    vector<Gene*> genes;
};

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

bool
is_compressed(const string &path){
    if (path.size() < 4){
        return false; // checking for .gz extension
    }
    uint idx = path.size() - 1;
    return path[idx - 2] == '.' && path[idx - 1] == 'g' && path[idx] == 'z'; 
}

void 
open_in_filestream(bool gzipped, gzFile &gz_fh, ifstream &fh, const string &infile){

    bool bad;
    if (gzipped){
        gz_fh = gzopen(infile.c_str(), "rb");
        bad   = gz_fh == NULL;
    }
    else {
        fh.open(infile);
        bad = !fh.is_open();
    }

    if (bad){
        cerr << "Error: could not open " << infile << '\n';
        exit(1);    
    }
}

void 
close_in_filestream(bool gzipped, gzFile &gz_fh, ifstream &fh){
    if (gzipped){
        gzclose(gz_fh);
    }
    else {
        fh.close();
    }
}

string
get_gzline(gzFile fh, bool &eof){
    //
    // construct a string that reaches the '\n' character
    //
    string line;
    const int buff_size = 8192;
    char buffer[buff_size];
    bool chars_read = false; // were characters read

    while (true){
        char *read_chars = gzgets(fh, buffer, buff_size);

        if (read_chars == NULL){
            break; // reach the end of the file stream
        }
        chars_read = true;
        line      += buffer;
        if (!line.empty() && line.back() == '\n'){
            break;
        }
    }

    eof = !chars_read; // will be true if no characters read
    return line;
}

int
parse_tabular(string &line, vector<string> &parts){

    //
    // parse a '\t' delimited line
    //

    int start  = 0, end = 0;

    //
    // start from an empty vector
    //
    parts.clear();

    while (end < line.size()){
        if (line[end] == '\t'){
            parts.emplace_back(line.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }

    if (start < line.size()){
        parts.emplace_back(line.substr(start));
    }

    return 0;
}

int
get_popmap(const string &infile, vector<string> &target){
    //
    // load the sample names for the target population
    // of interest
    //
    std::ifstream fh(infile);

    if (!fh.is_open()) {
        cerr << "Error: Could not open the file " << infile << '\n';
        exit(1);
    }

    string line;
    
    while (getline(fh, line)) {
        target.push_back(line);
    }

    fh.close();

    cerr << "Loaded " << target.size() << " samples for target population\n";

    return 0;
}

void
parse_attributes(string &attributes, vector<Attribute> &atrbVec){

    //
    // get the gene_id from a ';' delimited string
    // if we only want the gene_id
    //

    if (attributes.empty()){
        return;
    }

    size_t start = 0, next = string::npos, length = attributes.size();
    string part;

    atrbVec.clear(); // ensure new entries

    // iterate over a ';' delimited string
    while (start <= length){
        next = attributes.find(';', start);
        part = next == string::npos ? attributes.substr(start) : attributes.substr(start, next - start);
        
        // strip whitespace
        uint i = 0;
        while (i < part.size() && part[i] == ' '){
            i++;
        }

       part = part.substr(i);

       if (!part.empty()){
            size_t idx = part.find(' '); // find if we have a key value pair
            if (idx != string::npos){
                string key   = part.substr(0, idx);
                string value = part.substr(idx + 1);
                // remove qoutes
                if (value.size() >= 2 && value[0] == '"' && value.back() == '"'){
                    value = value.substr(1, value.size() - 2);
                }
                atrbVec.push_back({key, value});
            }
       }
        // we reached the end
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }
}

int
parse_annotation(const string &ann, unordered_map<string, vector<Gene*>> &genome, 
    const unordered_map<string, vector<SNP>> &markers){
    //
    // collect only genes on chromosomes with markers
    //

    bool gzipped = is_compressed(ann);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, ann);

    //
    // we will create a mapping
    // for gene_id -> Gene & at
    // the end, move them into 
    // the genome container
    //
    // chr -> gene_id -> Gene
    //
    unordered_map<string, unordered_map<string, Gene*>> gene_map;

    //
    // create the objects we need to store info
    //
    vector<string> parts;
    vector<Attribute> atrbVec;
    Gene *g;
    string line, chrom, gene_id, gene_name;
    uint start, end;
    bool eof = false;

    while (true){
        if (gzipped){
            line = get_gzline(gz_fh, eof);
            if (eof){
                break; // end of parsing
            } 
        }
        else {
            if (!getline(txt_fh, line)){
                break; // end of parsing
            }
        }

        if (eof){
            break; // end of file
        }
        if (line.empty() || line[0] == '#'){
            continue; // skip comment & empty lines
        }
        parse_tabular(line, parts);
        if (parts.size() < 9){
            cerr << "Malformed line in " << ann << '\n' << line;
            exit(1);
        }

        chrom = parts[0];

        if (markers.count(chrom) == 0){
            continue; // no markers on this chrom
        }

        if (parts[2] == "gene" || parts[2] == "exon"){

            if (parts[8].back() == '\n'){
                parts[8].pop_back(); // strip new line char
            }
            // grab the gene_id for this gene
            parse_attributes(parts[8], atrbVec);
            gene_id   = "";
            gene_name = "";
            for (uint i = 0; i < atrbVec.size(); i++){
                if (atrbVec[i].key == "gene_id"){
                    gene_id = atrbVec[i].val;
                }
                else if (atrbVec[i].key == "gene_name"){
                    gene_name = atrbVec[i].val;
                }
            }
            if (gene_id.empty()){
                cerr << "Unable to get gene_id for the following record:\n" << line;
                exit(1);
            }
            chrom = parts[0];
            start = (uint)stoi(parts[3]);
            end   = (uint)stoi(parts[4]);
            if (parts[2] == "gene"){
                g = new Gene(gene_id, gene_name, start, end);
                gene_map[chrom][gene_id] = g;
            }
            else if (gene_map[chrom].find(gene_id) != gene_map[chrom].end()) {
                Exon exon{start, end};
                gene_map[chrom][gene_id]->add_exon(exon);
            }
            else {
                cerr << "Malformed annotations. Exon came before gene entry. "
                     << "Offending line:\n" << line;
                exit(1);
            }
        } // end of exon parsing
    } // end of parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (gene_map.empty()){
        cerr << "No genes were found on marker-containing chromosomes in " << ann << '\n';
        exit(1);
    }
    // move the genes into the genome map

    for (auto itr = gene_map.begin(); itr != gene_map.end(); itr++){
        chrom                               = itr->first;
        unordered_map<string, Gene*> &genes = itr->second;
        vector<Gene*> &chrom_genes          = genome[chrom];

        chrom_genes.reserve(genes.size()); // reserve enough space
        for (auto jtr = genes.begin(); jtr != genes.end(); jtr++){
            jtr->second->resolve_exons(); // deduplicate & merge overlapping exons
            chrom_genes.push_back(jtr->second);
        }

        // now sort for downstream binary search
        sort(chrom_genes.begin(), chrom_genes.end(),[]
            (const Gene *geneA, const Gene *geneB){
                if (geneA->start == geneB->start){
                    return geneA->end < geneB->end;
                }
                return geneA->start < geneB->start;
            }  
        );
    }
    return 0;
}

void
get_overlapping_genes(const int left, const int site, vector<Gene*> *genes, vector<Gene *> &overlapped_genes){

    //
    // populate the candidates vector with overlapping genes
    // at this site
    //
    int right = left + 1;

    Gene *gene;

    while (right < genes->size()){
        gene = (*genes)[right];
        if (gene->start <= site && gene->end >= site){
            overlapped_genes.push_back(gene);
        }
        else if (gene->start > site){
            break;
        }
        right++;
    }
}

int
add_alleles(string &geno, vector<string> &alleles){
    
    // return 1 if this is a heterozygous
    // genotype, else 0

    if (geno.empty()){
        cerr << "Error: Encountered a not a valid genotype call\n";
        exit(1);
    }

    string a1, a2;
    uint start = 0, end = 0, length = geno.size();

    while (end < length){
        if (geno[end] == '/' || geno[end] == '|'){
            if (a1.empty()){
                a1 = geno.substr(start, end - start);
            }
            start = end + 1;
        }
        end++;
    }

    a2 = geno.substr(start);

    // add the alleles to the population's allele pool
    if (find(alleles.begin(), alleles.end(), a1) == alleles.end()){
        alleles.push_back(a1);
    }
    if (find(alleles.begin(), alleles.end(), a2) == alleles.end()){
        alleles.push_back(a2);
    }

    return 0;
}

int
get_markers(const string &vcf, vector<string> &target, 
            unordered_map<string, vector<SNP>> &markers){

    //
    // this function will be the main work horse for this
    // code
    //

    bool gzipped = is_compressed(vcf);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, vcf);
    
    //
    // sample -> number of total & hetero sites for this sample
    //
    
    vector<string> parts;
    vector<bool> in_target;
    string   line, sample, geno;
    bool     eof = false, found;

    uint total = 0, candidates = 0;

    while (true){
        if (gzipped){
            line = get_gzline(gz_fh, eof);
            if (eof){
                break; // end of parsing
            }
            
        }
        else {
            if (!getline(txt_fh, line)){
                break; // end of parsing
            }
        }

        if (line[0] == '#'){
            //
            // we may need to collect the sample indices
            // check for column starting with #CHROM 
            //
            if (line.size() > 6 && line.substr(0, 6) == "#CHROM"){
                if (line.back() == '\n'){
                    line.pop_back(); // we need to parse this file
                } 
                parse_tabular(line, parts);
                if (parts.size() < 10){
                    // there is no genotype info here
                    cerr << "No sample information found in " << vcf << '\n';
                    exit(1);
                }
                // now parse through
                for (uint i = 9; i < parts.size(); i++){
                    sample = parts[i];
                    found  = false;
                    for (uint j = 0; j < target.size(); j++){
                        if (sample == target[j]){
                            found = true;
                            break;
                        }
                    }
                    in_target.push_back(found);
                }
            }
            continue;
        }

        if (in_target.empty()){
            cerr << "Error: Could not find #CHROM header containg sample information\n";
            exit(1);
        }

        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }

        total++;
        parse_tabular(line, parts);

        // first, check that the site is homozygous
        // for the target population
        vector<string> target_alleles;
        for (uint i = 9; i < parts.size(); i++){
            uint sample_idx = i - 9;
            if (!in_target[sample_idx]){
                continue;
            }
            geno = parts[i];
            geno = geno.substr(0, geno.find(':')); 
            if (geno.find('.') != string::npos){ // missing data
                continue;
            }
            add_alleles(geno, target_alleles);
            if (target_alleles.size() > 1){ // heterozgyous
                break;
            }
        }

        if (target_alleles.size() != 1){
            continue; // heterzygous or no information
        }

        // get genotypes of other samples
        vector<string> other_alleles;
        const char t_allele = target_alleles[0].at(0);
        bool valid = true;
        for (uint i = 9; i < parts.size(); i++){
            uint sample_idx = i - 9;
            if (in_target[sample_idx]){
                continue;
            }
            geno = parts[i];
            geno = geno.substr(0, geno.find(':'));

            // skip missing data & verify that it is not homozygous
            // for the target allele
            uint count   = 0;
            bool missing = false;
            for (uint j = 0; j < geno.size(); j++){
                if (geno[j] == '.'){
                    missing = true;
                    break;
                }
                else if (geno[j] == t_allele){
                    count++;
                }
            }
            if (missing){
                continue;
            }
            else if (count > 1){
                valid = false;
                break;
            }
            add_alleles(geno, other_alleles);
        }

        if (!valid || other_alleles.empty()){
            continue;
        }

        SNP snp{(uint)stoi(parts[1])};
        markers[parts[0]].push_back(snp);
        candidates++;


    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (total == 0){
        cerr << "Found no records in this file\n";
        exit(1);
    }

    double p = ((double)candidates/ (double)total) * 100.0;
    cerr << std::fixed << setprecision(2);
    cerr << "Scanned " << total << " sites and found " <<  candidates << " (" << p << "%) candidates\n";

    return 0;
}

int
overlap_genes(unordered_map<string, vector<Gene*>> &genome,
                   unordered_map<string, vector<SNP>> &markers){
    // add each gene to any SNPs that falls within its boundaries

    uint total = 0;
    for (auto itr = markers.begin(); itr != markers.end(); itr++){
        
        // nothing to overlap here
        if (genome.find(itr->first) == genome.end()){
            continue;
        }
        string chrom         = itr->first;
        vector<SNP> &snps    = itr->second;
        vector<Gene*> *genes = &genome[chrom];

        // genes are already sorted, but we need to keep
        // track when a gene ends
        uint furthest = 0, num_genes = genes->size();
        vector<uint> gene_ends(num_genes, 0);
        for (uint i = 0; i < num_genes; i++){
            if (furthest < (*genes)[i]->end){
                furthest = (*genes)[i]->end;
            }
            gene_ends[i] = furthest;
        }

        for (uint i = 0; i < snps.size(); i++){
            SNP &snp = snps[i];

            // use a binary search to find the leftmost
            // gene this snp can overlap
            Gene* gene;
            uint pos = snp.pos;
            uint left = 0, right = num_genes, mid;

            while (left < right){
                mid  = left + (right - left) / 2;
                gene = (*genes)[mid];
                if (gene_ends[mid] >= pos){
                    right = mid;
                }
                else {
                    left = mid + 1;
                }
            }

            bool overlaps = false;
            if (left < num_genes){
                gene     = (*genes)[left];
                overlaps = gene->start <= pos;
            }
            if (overlaps){
                total++;
                vector<Gene *> overlapped_genes = {gene};
                get_overlapping_genes(left, pos, genes, overlapped_genes);
                snp.genes = overlapped_genes;
            }

        } // end of snps loop
    } // end of chrom loop

    cerr << "Found " << total << " snps that overlap genes\n";

    return 0;

}

int
clear_genes(unordered_map<string, vector<Gene*>> &genome){
    // helper function to clear the gene objects

    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene*> &genes = itr->second;
        for (uint i = 0; i < genes.size(); i++){
            delete genes[i];
        }
    }
    return 0;
}

int
write_output(unordered_map<string, vector<SNP>> &markers){
    //
    // write each SNP and the genes it overlaps
    //

    ofstream fh("homozygous_sites.tsv");

    if (!fh.is_open()){
        cerr << "Error: Unable to open an output file in this directory\n";
        exit(1);
    }

    fh << "#Chrom\tBP\tGenes\n";

    for (auto itr = markers.begin(); itr != markers.end(); itr++){
        string chrom      = itr->first;
        vector<SNP> &snps = itr->second;

        for (uint i = 0; i < snps.size(); i++){
            SNP &snp              = snps[i];
            vector<Gene *> &genes = snp.genes;
            fh << chrom << '\t' << snp.pos << '\t';

            // if it's empty, this will not run
            for (uint j = 0; j < genes.size(); j++){
                Gene* gene  = genes[j];
                SnpType t   = gene->classify(snp.pos);
                string type = t == SnpType::Exonic ? " exonic" : " intronic";
                string name = gene->name.empty() ? gene->id : gene->name;
                if (j > 0){
                    fh << ", ";
                }
                fh << name << type;
            } // end of genes loop
            fh << '\n';
        } // end of snps loop
    } // end of markers loop

    fh.close();

    return 0;
}

void
help(){
    cerr << "Usage: ./find_homo_sites -v vcf_file -p pop_map -a ann.gtf\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, popmap;
    
    // expect at least 3 arguments
    if (argc < 7){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = string(argv[i + 1]);
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = string(argv[i + 1]);
        }
        else if (arg == "-p" && i + 1 < argc){
            popmap = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (vcf.empty() && ann.empty() && popmap.empty()){
        help();
    }

    if (!file_exists(vcf)){
        cerr << "Unable to find " << vcf << '\n';
        exit(1);
    }
    if (!file_exists(ann)){
        cerr << "Unable to find " << ann << '\n';
        exit(1);
    }
    if (!file_exists(popmap)){
        cerr << "Unable to find " << popmap << '\n';
        exit(1);
    }

    // load the target population
    vector<string> pop;
    get_popmap(popmap, pop);

    // load candidate markers
    unordered_map<string, vector<SNP>> markers;
    get_markers(vcf, pop, markers);

    // grab genes on focal chromosomes
    unordered_map<string, vector<Gene*>> genome;
    parse_annotation(ann, genome, markers);

    // overlap the two dataset
    overlap_genes(genome, markers);

    // write the output
    write_output(markers);

    // clean up
    clear_genes(genome);

    return 0;
}