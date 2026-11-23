#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <zlib.h>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <cmath>
#include <random>

using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::sort;
using std::find;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::setprecision;
using std::stod;
using std::stoi;

//
// this code is specifically designed. for finding heterozygous
// sites for the population of interest (bear paw). If given an
// annotation file, only exonic sides will be searched. If given a
// list of genes, then only those genes within the annotation will be
// searched
//

typedef unsigned int uint;

//
// LD analysis was adapted from Kuster et al. (2025)
//  Kuster, S. A., Schumer, M., Havird, J., & Sloan, D. B. (2026). 
//  Sensitivity of genome-wide tests for mitonuclear genetic incompatibilities. 
//  Philosophical transactions of the Royal Society of London. Series B, Biological sciences, 381(1947), 20250084. 
//
//
// mitonuclear linkage disequilibrium (LD) within bear paw
//
// A = the nuclear ENA allele (defined by the pure populations)
// B = the ENA mitotype (from the popmap)
// every nuclear allele copy is paired with its fish's mitotype,
// so no phasing is needed
//
//   D  = pAB - pA * pB
//   D' = D / Dmax
//   r2 = D^2 / (pA (1 - pA) pB (1 - pB))
//
// D > 0 -> matched (parental) combinations are more common
// D < 0 -> mismatched combinations are more common
//

const double MAX_EXACT_PERMS = 100000; // enumerate every relabeling up to this many
const uint   N_RANDOM_PERMS  = 10000;  // otherwise use random relabelings
const double PERM_TOL        = 1e-12;  // tolerance when comparing |D| values

struct MitoLD {
    uint   n_ena     = 0;  // allele copies in bear paw ENA fish
    uint   n_tnp     = 0;  // allele copies in bear paw TNP fish
    double q_ena     = 0;  // ENA allele freq in bear paw ENA fish
    double q_tnp     = 0;  // ENA allele freq in bear paw TNP fish
    double p_b       = 0;  // fraction of allele copies from ENA mitotype fish
    double d         = 0;
    double d_prime   = 0;
    double r2        = 0;
    double p_perm    = 1;  // two sided permutation p-value for D
    uint   n_perm    = 0;  // number of relabelings used
    bool   exact     = true;
    bool   valid     = false; // can LD be computed at all
    bool   dp_valid  = false; // is D' defined
    bool   r2_valid  = false; // is r2 defined

    // heterozygosity in all diploid bear paw calls
    uint   n_dip     = 0;
    double ho        = 0;
    double he        = 0;
    double fis       = 0;
    bool   fis_valid = false;
};


enum Sex {
    Male   = 0,
    Female = 1
};

enum Strand {
    Unknown = 0,
    Forward = 1,
    Reverse = 2
};

enum mtType {
    ENA = 0,
    TNP = 1
};

enum RegionType {
    NoRegion = 0,
    Exonic   = 1,
    Promoter = 2
};

enum Classification {
    Uinform   = 0, // un-informative
    Diff      = 1,
    Private   = 2,
    BpHetSame = 3, // pure pops fixed for the same allele, bear paw het
    BpHetDiff = 4, // pure pops fixed for different alleles, bear paw het
};
struct Exon {
    uint start;
    uint end;
};
struct Attribute {
    string key;
    string val;
};

struct Sample {
    mtType pop;
    Sex    sex;
    bool   bear_paw;
};

class Gene {
private:
    bool _coding = false;

public:
    string       id;
    string       name;
    uint         start;
    uint         end;
    Strand       strand;
    vector<Exon> exons;

    //
    // members for searching for sites possibly
    // upstream in the promoter region
    //
    bool has_promoter = false;
    uint prom_start   = 0;
    uint prom_end     = 0;
    uint search_start = 0;
    uint search_end   = 0;

    Gene (string id_, string name, uint start, uint end, Strand strand){
        this->id           = id_;
        this->name         = name;
        this->start        = start;
        this->end          = end;
        this->strand       = strand;
        this->search_start = start;
        this->search_end   = end;
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

    bool is_exonic(const uint pos){
        //
        // exons are ordered,
        // so just see if there is an overlap
        //

        for (uint i = 0; i < this->exons.size(); i++){
            const Exon &e = this->exons[i];
            if (e.start <= pos && pos <= e.end){
                return true;
            }
        }
        // no overlap
        return false;
    }

    bool in_promoter(const uint pos){
        //
        // is this site in the promoter search space
        //
        if (this->has_promoter){
            return this->prom_start <= pos && pos <= this->prom_end;
        }
        return false;
    }

    void set_promoter(const uint upstream){

        //
        // the promoter is upstream of the first
        // exon accounting for strand orientation
        //

        if (upstream == 0 || this->strand == Strand::Unknown){
            return; // nothing to do
        }

        if (this->strand == Strand::Forward){
            if (this->start <= 1){
                return; // edge case
            }
            if (this->start > upstream){ // avoid negative positions
                this->prom_start = this->start - upstream;
            }
            else {
                this->prom_start = 1;
            }
            this->prom_end     = this->start - 1;
            this->search_start = this->prom_start;
        }
        else if (this->strand == Strand::Reverse){
            this->prom_start = this->end + 1;
            this->prom_end   = this->end + upstream;
            this->search_end = this->prom_end;
        }
        else {
            return; // unknown strand, shouldn't trigger
        }

        this->has_promoter = true;
    }

    void set_as_coding(void){
        this->_coding = true;
    }

    bool is_coding(void) const {
        return this->_coding;
    }
 
};

//
// global constants for this analysis
//
const string BP_PREFIX = "Bear_Paw"; // sample ids starting with this are bear paw
const string MT        = "MT";       // haploid for all
const string Y         = "Y";        // haploid for males, missing in females
const string X         = "XIX";      // haploid for males 2.5 Mb and onwards
const uint   PAR_END   = 2500000;    // end of the pseudoautosomal region on X

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

bool
read_line(bool gzipped, gzFile &gz_fh, ifstream &fh, string &line){
    //
    // return false at the end of the file
    // \n chars are removed
    //

    bool eof;
    if (gzipped){
        line = get_gzline(gz_fh, eof);
        if (eof){
            return false;
        }
    }
    else if (!getline(fh, line)){
        return false;
    }

    // we didn't return, so strip \n
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
        line.pop_back();
    }

    return true;
}


int
parse_popmap(const string &pop_file, unordered_map<string, Sample*> &pop_map, vector<string> &pops){

    //
    // read in the tsv file that contains the population map
    //
    //

    ifstream fh(pop_file);

    if (!fh.is_open()){
        cerr << "Error: Unable to open " << pop_file << '\n';
        exit(EXIT_FAILURE);
    }

    string line, sample, pop;
    vector<string> parts;
    uint males = 0, females = 0;
    bool bear_paw;
    char sex;

    while (getline(fh, line)){

        while (!line.empty() && line.back() == '\r'){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }

        parse_tabular(line, parts);

        if (parts.size() != 3){
            cerr << "Error: Invalid format in mapping file. Expected 3 columns."
                 << " Each row should be [sample][pop][sex]."
                 << " Offending line:\n" << line;

            exit(EXIT_FAILURE);
        }
        sample = parts[0];
        pop    = parts[1];
        sex    = parts[2][0];

        if (sex != 'M' && sex != 'F'){
            cerr << "Error: Found an invalid Sex entry. Expected M or F in column 3 in pop map."
                 << " Offending line:\n" << line;
            exit(EXIT_FAILURE);
        }
        bear_paw = sample.compare(0, BP_PREFIX.size(), BP_PREFIX) == 0;
        mtType m = pop == "ENA" ? mtType::ENA : mtType::TNP;
        Sex    s = sex == 'M' ? Sex::Male : Sex::Female;
        pop_map[sample] = new Sample{m, s, bear_paw};

        if (find(pops.begin(), pops.end(), pop) == pops.end()){
            pops.push_back(pop);
        }
        if (sex == 'M'){
            males++;
        }
        else {
            females++;
        }
    }

    fh.close();

    if (pop_map.empty()){
        cerr << "No samples were found in the population map\n";
        exit(EXIT_FAILURE);
    }
    if (pops.size() != 2){
        cerr << "Error: Expected exactly 2 populations in the pop map. Found " << pops.size() << '\n';
        exit(EXIT_FAILURE);
    }

    cerr << "Loaded " << pop_map.size() << " samples across " << pops.size() << " populations\n";
    cerr << males << " were found to be males while the remaining (" << females << ") were found to be female.\n";

    // ensure ENA comes before TNP
    sort(pops.begin(), pops.end());

    return 0;

};

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
                 const bool only_targets, unordered_set<string> &target_genes, uint upstream){
    //
    // collect genes in this specific file and optionally only
    // those in the targets container
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
    bool keep_reading;

    while (true){
        keep_reading = read_line(gzipped, gz_fh, txt_fh, line);

        if (keep_reading == false){
            break;
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

        if (parts[2] == "gene" || parts[2] == "exon" || parts[2] == "CDS"){

            if (parts[8].back() == '\n'){
                parts[8].pop_back(); // strip new line char
            }
            // grab the gene_id for this gene
            parse_attributes(parts[8], atrbVec);
            gene_id.clear();
            gene_name.clear();
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
            else if (only_targets && target_genes.count(gene_id) == 0){
                continue; // skip this
            }
            chrom = parts[0];
            start = (uint)stoi(parts[3]);
            end   = (uint)stoi(parts[4]);
            if (parts[2] == "gene"){
                Strand strand;
                switch (parts[6][0]) {
                    case '+':
                        strand = Strand::Forward;
                        break;
                    case '-':
                        strand = Strand::Reverse;
                        break;
                    default:
                        strand = Strand::Unknown;
                        break;
                    }
                g = new Gene(gene_id, gene_name, start, end, strand);
                gene_map[chrom][gene_id] = g;
            }
            else if (gene_map[chrom].find(gene_id) != gene_map[chrom].end()) {
                if (parts[2] == "exon"){
                    Exon exon{start, end};
                    gene_map[chrom][gene_id]->add_exon(exon);
                }
                else {
                    gene_map[chrom][gene_id]->set_as_coding();
                }
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
        cerr << "No genes were found in " << ann << '\n';
        exit(1);
    }

    // move the genes into the genome map
    uint total = 0;
    for (auto itr = gene_map.begin(); itr != gene_map.end(); itr++){
        chrom                               = itr->first;
        unordered_map<string, Gene*> &genes = itr->second;
        vector<Gene*> &chrom_genes          = genome[chrom];

        chrom_genes.reserve(genes.size()); // reserve enough space
        for (auto jtr = genes.begin(); jtr != genes.end(); jtr++){

            if (!jtr->second->is_coding()){
                delete jtr->second; // free the memory
                continue;           // skip non-coding genes
            }
            jtr->second->resolve_exons(); // collapse & merge exons
            jtr->second->set_promoter(upstream);
            chrom_genes.push_back(jtr->second);
            total++;
        }

        // now sort for downstream binary search
        sort(chrom_genes.begin(), chrom_genes.end(),[]
            (const Gene *geneA, const Gene *geneB){
                if (geneA->search_start == geneB->search_start){
                    return geneA->search_end < geneB->search_end;
                }
                return geneA->search_start < geneB->search_start;
            }  
        );
    }

    cerr << "Loaded " << total << " genes\n";
    return 0;
}

int
load_target_genes(const string &genes_file, const bool only_targets,
                  unordered_set<string> &targets){

    // helper function to create a container holding the gene_ids of interest

    if (only_targets == false){
        return 1; // we are not using this feature
    }

    ifstream fh(genes_file);

    if (!fh.is_open()) {
        cerr << "Error: Could not open the file " << genes_file << '\n';
        exit(1);
    }

    string line;
    
    while (getline(fh, line)) {
        targets.insert(line);
    }

    fh.close();

    cerr << "Loaded " << targets.size() << " target genes\n";

    return 0;
}

int
clean_up(unordered_map<string, vector<Gene*>> &genome, unordered_map<string, Sample*> &pop_map){
    // just some code to clear up the memory for each pointer

    for (auto itr = genome.begin(); itr != genome.end(); itr++){

        vector<Gene*> &genes = itr->second;
        for (uint i = 0; i < genes.size(); i++){
            delete genes[i];
        }
    }

    for (auto itr = pop_map.begin(); itr != pop_map.end(); itr++){
        delete itr->second;
    }

    return 0;
}

Gene *
is_candidate_region(const uint pos, vector<Gene*>*genes, vector<uint> &gene_ends, RegionType &region){

    //
    // helper function to determine
    // whether this site lands on an exonic region, a candidate promoter region, 
    // or neither
    //

    region = RegionType::NoRegion;

    if (genes->empty()){
        return nullptr; // no genes to search through
    }

    const uint ngenes = genes->size();

    // binary search for a hit
    uint left = 0, mid, right = ngenes;

    while (left < right){
        mid = left + (right - left) / 2;
        if (gene_ends[mid] >= pos){
            right = mid;
        }
        else {
            left = mid + 1;
        }
    }

    Gene *gene;
    Gene *promoter_hit = nullptr;

    while (left < ngenes){
        gene = (*genes)[left];
        
        if (pos < gene->search_start){
            break; // this & all later genes begin after pos
        }
        if (pos <= gene->search_end){
            if (gene->start <= pos && pos <= gene->end && gene->is_exonic(pos)){
                region = RegionType::Exonic;
                return gene; // lands in an exon
            }
            else if (promoter_hit == nullptr && gene->in_promoter(pos)){
                promoter_hit = gene; // keep looking for an exon
            }
        }
        left++;
    }

    if (promoter_hit != nullptr){
        region = RegionType::Promoter;
    }

    return promoter_hit; // this may never change from nullptr
}

bool
is_biallelic(const string &ref, const string &alt){
    //
    // has there been 2 alleles reported at this site
    //

    if (ref.size() != 1 || alt.size() != 1){
        return false;
    }
    if (alt == "*" || alt == "."){
        return false; // spans deletion or no alt
    }
    return true;
}

int
expected_ploidy(const string &chrom, const uint pos, const Sex sex){
    //
    // return the number of expected alleles at this site
    //

    if (chrom == MT){
        return 1; // all samples should be haploid
    }
    if (chrom == Y){
        return sex == Sex::Male ? 1 : 0;
    }
    if (chrom == X && pos > PAR_END && sex == Sex::Male){
        return 1; // only applies to males
    }

    return 2;
}

int
count_allele_copies(const string &field, const uint gt_idx, const uint ploidy, const int allele){

    //
    // count the copies of allele (0 == ref, 1 = alt) for this genotype
    // return -1 when the call is missing or does not match the expected
    // ploidy
    //

    if (ploidy == 0){
        return -1; // no call should be made here
    }

    //
    //
    // walk to the GT subfield
    //
    uint start = 0, idx = 0, length = field.size();

    while (idx < gt_idx && start < length){
        if (field[start] == ':'){
            idx++;
        }
        start++;
    }

    if (idx < gt_idx){
        cerr << "Encountered a truncated field\n";
        return -1; // should not happen (truncated sample field)
    }

    //
    // collect the alleles by constructing
    // an array of size 2 (assumes biallelic sites)
    //
    int alleles[2];
    uint ncalled = 0;

    for (uint i = start; i < length && field[i] != ':'; i++){
        char c = field[i];
        if (c == '/' || c == '|'){
            continue;
        }
        if (c != '0' && c != '1'){
            return -1; // missing ('.') or unexpected allele
        }
        if (ncalled == 2){
            cerr << "Encountered a sample with >2 alleles\n";
            return -1; // should not happen
        }
        alleles[ncalled++] = c - '0';
    } // end of i

    if (ncalled == 0){
        cerr << "Encountered a sample with no alleles but expected?\n";
        return -1; // 
    }

    if (ploidy == 1){
        if (ncalled == 2 && alleles[0] != alleles[1]){
            cerr << "Encountered a heterozygous call in a haploid region\n";
            return -1;
        }
        return alleles[0] == allele ? 1 : 0;
    }
    if (ncalled != 2){
        cerr << "Encountered a haploid call in a heterozygous region\n";
        return -1;
    }

    //
    // count the number of copies
    //
    int copies = 0;
    for (uint i = 0; i < 2; i++){
        if (alleles[i] == allele){
            copies++;
        }
    }

    return copies;

}

bool
passes_missing(const vector<Sample*> &samples, const vector<uint> &ploidy,
               const vector<int> &alt_copies, const double max_missing){

    //
    // determine whether the current record has enough samples (based on max_missing)
    // to consider this informative/useful
    //
    // for the Y chrom, females will not be penalized
    // 
    // bp = Bear Paw
    // expected[pop][bp] & missing[pop][bp]
    // pop: 0 = ENA, 1 = TNP, bp: 0 = "pure", 1 = bear paw
    //
    uint expected[2][2] = {{0, 0}, {0,0}};
    uint missing[2][2]  = {{0, 0}, {0,0}};

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];

        if (sample == nullptr || ploidy[i] == 0){
            continue;
        }
        uint bp = 0;
        if (sample->bear_paw == true){
            bp = 1;
        }
        expected[sample->pop][bp]++;

        if (alt_copies[i] < 0){
            missing[sample->pop][bp]++;
        }

    } // end of i

    // check if we passed the thresold
    for (uint i = 0; i < 2; i++){
        for (uint j = 0; j < 2; j++){
            if (expected[i][j] == 0){
                continue;
            }
            double frac = (double)missing[i][j] / expected[i][j];

            if (frac > max_missing){
                return false;
            }
        }
    }

    // we never passed the threshold
    return true;
}

void
tally_genotypes(const vector<string> &parts, const vector<Sample*> &samples, const int gt_idx,
                const string &chrom, const uint pos, vector<uint> &ploidy, vector<int> &alt_copies,
                uint pure_alt[2], uint pure_n[2]){

    //
    // fill the pure matrices for every sample
    // the idea is that one population should have counts for one allele
    // while the other population should have counts for the alternative allele
    // index 0 = ENA, 1 = TNP
    //
    
    // set the counts to 0
    for (uint i = 0; i < 2; i++){
        pure_alt[i] = 0;
        pure_n[i]   = 0;
    }

    // now loop through all the samples
    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];

        if (sample == nullptr){
            continue; // not in pop map
        }
        ploidy[i]     = expected_ploidy(chrom, pos, sample->sex);
        alt_copies[i] = count_allele_copies(parts[i + 9], gt_idx, ploidy[i], 1);

        // negative alt_copies[i] (i.e., -1) means the site is missing
        if (sample->bear_paw == true || alt_copies[i] < 0){
            continue;
        }

        //
        // pure_n is counting the number of alleles for each pop
        // so a diploid individual will contribute 2 alleles
        // needed since males & females can differ 
        //

        pure_alt[sample->pop] += alt_copies[i];
        pure_n[sample->pop]   += ploidy[i];
    }
}

void
tally_bear_paw(const vector<Sample *> &samples, const vector<uint> &ploidy,
               const vector<int> &alt_copies, const int a_allele, uint bp[2][4]){
    
    //
    // count the genotypes for bear paw samples in terms of
    // ENA [0] or TNP [1] alleles
    //
    // bp[group]: [0] homo ENA, [1] het, [2] homo TNP, [3] missing
    // haploid calls will be treated as homo
    //

    // prefill with 0's
    for (uint i = 0; i < 2; i++){
        for (int j = 0; j < 4; j++){
            bp[i][j] = 0;
        }
    }

    // now go through each bear paw sample
    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i]; // ploidy == 0 can happen for females at Y
        if (sample == nullptr || !sample->bear_paw || ploidy[i] == 0){
            continue;
        }

        if (alt_copies[i] < 0){
            // -1 count, so possibly missing
            bp[sample->pop][3]++;
            continue;
        }

        // count ENA copies
        int ena_copies;
        if (a_allele == 1){
            ena_copies = alt_copies[i];
        }
        else {
            ena_copies = (int)ploidy[i] - alt_copies[i];
        }

        if (ena_copies == int(ploidy[i])){
            bp[sample->pop][0]++; // homo ENA
        }
        else if (ena_copies == 0){
            bp[sample->pop][2]++; // homo TNP
        }
        else {
            bp[sample->pop][1]++; // het
        }
    } // end of i
}

int
classify_site(const string &chrom, const uint pure_alt[2], const uint pure_n[2],
              const double min_freq, Classification &c){

    //
    // classify whether the site is differentiating or there is an informative
    // private allele. 
    // return 0  -> ref belongs to ENA
    // return 1  -> alt belongs to ENA
    // return -1 -> non-informative
    //
    
    if (pure_n[0] == 0 || pure_n[1] == 0){
        return -1; // no calls ine one of the "pure" populations
    }

    double alt_freq_a      = (double)pure_alt[0] / pure_n[0];
    double alt_freq_b      = (double)pure_alt[1] / pure_n[1];
    const bool haploid     = chrom == MT || chrom == Y;
    const double site_freq = haploid ? 1.0 : min_freq; // haploid sites must be fixed

    if (alt_freq_a >= site_freq && alt_freq_b <= 1.0 - site_freq){
        c = Classification::Diff;
        return 1;
    }
    else if (alt_freq_b >= site_freq && alt_freq_a <= 1.0 - site_freq){
        c = Classification::Diff;
        return 0;
    }
    else if (haploid){
        return -1; // haploid call should have been fixed
    }

    // if none returned, check for a private allele
    if (pure_alt[0] > 0 && pure_alt[1] == 0){
        c = Classification::Private;
        return 1; // alt only seen in ENA
    }
    else if (pure_alt[1] > 0 && pure_alt[0] == 0){
        c = Classification::Private;
        return 0; // only seen in TNP
    }

    // checking if ref allele is private
    uint ref_a = pure_n[0] - pure_alt[0]; // num of total alleles - num of alt alleles
    uint ref_b = pure_n[1] - pure_alt[1];

    if (ref_a > 0 && ref_b == 0){
        c = Classification::Private;
        return 0; // ref only seen in ENA
    }
    else if (ref_b > 0 && ref_a == 0){
        c = Classification::Private;
        return 1; // ref only seen in TNP
    }

    return -1;
}

int
classify_fixed_pure(const vector<Sample*> &samples, const vector<uint> &ploidy,
                    const vector<int> &alt_copies, Classification &c){

    //
    // helper function to determine if the bear paw samples are hetero while
    // the other pure populations remain homozygous regardless of allele
    //
    // fixed[pop]: -1 not seen yet, 0 = ref, 1 = alt
    //

    int fixed[2] = {-1, -1};

    // tally pure pops only first
    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];

        if (sample == nullptr || sample->bear_paw || ploidy[i] == 0 || alt_copies[i] == -1){
            continue;
        }

        int allele;
        if (alt_copies[i] == 0){
            allele = 0; // homo ref
        }
        else if (alt_copies[i] == (int)ploidy[i]){
            allele = 1; // homo alt
        }
        else {
            return -1; // this individual is hetero
        }

        if (fixed[sample->pop] == -1){
            fixed[sample->pop] = allele; // assign
        }
        else if (fixed[sample->pop] != allele){
            return -1; // both homozygotes within one population
        }
    }

    if (fixed[0] == -1 || fixed[1] == -1){
        return -1; // no calls
    }

    if (fixed[0] == fixed[1]){
        c = Classification::BpHetSame;
    }
    else {
        c = Classification::BpHetDiff; // would be ideal
    }

    return fixed[0];
}

bool
bear_paw_all_het(const vector<Sample*> &samples, const vector<uint> &ploidy,
                 const vector<int> &alt_copies, const double min_frac){
    //
    // helper function to determine if all the bear paw samples are
    // hetero at this site 
    //

    uint called[2] = {0, 0};
    uint het[2]    = {0, 0};

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || !sample->bear_paw || ploidy[i] != 2 || alt_copies[i] == -1){
            continue; // only diploid bear paw samples, haploid, or missing
        }

        called[sample->pop]++;
        if (alt_copies[i] == 1){
            het[sample->pop]++;
        }
    }

    for (uint i = 0; i < 2; i++){
        if (called[i] == 0){
            return false; // can't evaluate
        }
        double frac = (double)het[i] / called[i];
        if (frac < min_frac){
            return false;
        }
    }

    return true;
}

int
find_gt_index(const string &format){

    //
    // find which ':' delimited subfield
    // holds the genotype (GT) info
    //

    uint start = 0, idx = 0, length = format.size();

    for (uint i = 0; i <= length; i++){
        if (i == length || format[i] == ':'){
            if (format.compare(start, i - start, "GT") == 0){
                return idx;
            }
            idx++;
            start = i + 1;
        }
    }

    // no GT info found
    return -1;
}

bool
shares_homozygotes(const vector<Sample*> &samples, const vector<uint> &ploidy,
                   const vector<int> &alt_copies){

    //
    // helper function to enfoce that the TNP and ENA populations have a homozygous
    // individual for the same allele. For example, both TNP and ENA can have an A/A individual
    //

    uint hom_ref[2] = {0, 0};
    uint hom_alt[2] = {0, 0};

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];

        if (sample == nullptr || sample->bear_paw || ploidy[i] == 0){
            continue;
        }

        if (alt_copies[i] == -1){
            continue; // missing data
        }

        if (alt_copies[i] == 0){
            hom_ref[sample->pop]++;
        }
        else if (alt_copies[i] == (int)ploidy[i]){
            hom_alt[sample->pop]++;
        }
    } // end of i

    if (hom_ref[0] > 0 && hom_ref[1] > 0){
        return true;
    }
    else if (hom_alt[0] > 0 && hom_alt[1] > 0){
        return true;
    }

    return false;
}

bool
is_informative(const uint bp[2][4], const double min_frac, const Classification c,
               const int ena_allele, const uint pure_alt[2], const uint pure_n[2]){

    //
    // helper function to determine if this site is informative
    // Min fract of Bear Paw ENA must have >= 1 ENA allele
    // Min fract of Bear Paw TNP must have >= 1 TNP allele
    //

    uint called_ena = bp[0][0] + bp[0][1] + bp[0][2];
    uint called_tnp = bp[1][0] + bp[1][1] + bp[1][2];

    if (called_ena == 0 || called_tnp == 0){
        return false; // can't evaluate
    }

    // hom or het for their type's allele counts as carrying it
    uint ena_carriers = bp[0][0] + bp[0][1];
    uint tnp_carriers = bp[1][2] + bp[1][1];

    //
    // be stricter on private alleles. If the genotypes
    // between the two "pure" populations can be found, 
    // we will not consider the site anymore
    //

    if (c == Classification::Private){
         //
        // find which population holds the private allele
        // the other population has 0 copies of it
        //
        uint ena_copies_in_tnp; // copies of the ENA allele in pure TNP
        uint tnp_copies_in_ena; // copies of the TNP allele in pure ENA

        if (ena_allele == 1){ // alt is ENA's allele
            ena_copies_in_tnp = pure_alt[1];
            tnp_copies_in_ena = pure_n[0] - pure_alt[0];
        }
        else { // ref is ENA's allele
            ena_copies_in_tnp = pure_n[1] - pure_alt[1];
            tnp_copies_in_ena = pure_alt[0];
        }

        if (tnp_copies_in_ena == 0){
            // private allele is in TNP, so BP ENA must be hom ENA
            ena_carriers = bp[0][0];
        }
        else if (ena_copies_in_tnp == 0){
            // private allele is in ENA, so BP TNP must be hom TNP
            tnp_carriers = bp[1][2];
        }
    }

    double frac_ena = (double)ena_carriers / called_ena;
    double frac_tnp = (double)tnp_carriers / called_tnp;

    return (frac_ena >= min_frac) && (frac_tnp >= min_frac);
}

double
ld_d(const int group_a, const int group_n, const int total_a, const int total_n){

    //
    // calculate D (D = pAB - pA * pB)
    // group_a: ENA allele copies in the fish labeled ENA mitotype
    // group_n: all allele copies in the fish labeled ENA mitotype
    //

    double p_a  = (double)total_a / total_n;
    double p_b  = (double)group_n / total_n;
    double p_ab = (double)group_a / total_n;

    return p_ab - p_a * p_b;
}

double
n_choose_k(const uint n, const uint k){

    //
    // number of ways to choose k of n, as a double to avoid overflow
    //
    double result = 1;
    for (uint i = 1; i <= k; i++){
        result = result * (n - k + i) / i;
    }
    return result;
}

void
mitonuclear_ld(const vector<Sample*> &samples, const vector<uint> &ploidy,
               const vector<int> &alt_copies, const int ena_allele, MitoLD &ld){

    //
    // compute the LD between the nuclear ENA allele & the ENA mitotype
    // in the bear paw samples, along with a permutation test where
    // the mitotype labels are shuffled among the called bear paw fish
    //

    ld = MitoLD(); // reset

    // collect the called bear paw fish
    vector<int>  a_copies; // ENA allele copies per fish
    vector<int>  n_copies; // allele copies per fish (ploidy)
    vector<bool> is_ena;   // ENA mitotype?
    uint het = 0, dip_alt = 0;

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || !sample->bear_paw || ploidy[i] == 0){
            continue;
        }
        if (alt_copies[i] < 0){
            continue; // missing
        }

        int ena_copies;
        if (ena_allele == 1){
            ena_copies = alt_copies[i];
        }
        else {
            ena_copies = (int)ploidy[i] - alt_copies[i];
        }

        a_copies.push_back(ena_copies);
        n_copies.push_back((int)ploidy[i]);

        if (sample->pop == mtType::ENA){
            is_ena.push_back(true);
            ld.n_ena += ploidy[i];
        }
        else {
            is_ena.push_back(false);
            ld.n_tnp += ploidy[i];
        }

        // heterozygosity only uses diploid calls
        if (ploidy[i] == 2){
            ld.n_dip++;
            dip_alt += alt_copies[i];
            if (alt_copies[i] == 1){
                het++;
            }
        }
    }

    // observed & expected heterozygosity
    if (ld.n_dip > 0){
        double p = (double)dip_alt / (2 * ld.n_dip);
        ld.ho    = (double)het / ld.n_dip;
        ld.he    = 2 * p * (1 - p);
        if (ld.he > 0){
            ld.fis       = 1 - ld.ho / ld.he;
            ld.fis_valid = true;
        }
    }

    if (ld.n_ena == 0 || ld.n_tnp == 0){
        return; // need both mitotypes
    }

    // observed sums
    const uint nfish = a_copies.size();
    int total_a = 0, total_n = 0, ena_a = 0, ena_k = 0;
    for (uint i = 0; i < nfish; i++){
        total_a += a_copies[i];
        total_n += n_copies[i];
        if (is_ena[i]){
            ena_a += a_copies[i];
            ena_k++;
        }
    }

    ld.valid = true;
    ld.q_ena = (double)ena_a / ld.n_ena;
    ld.q_tnp = (double)(total_a - ena_a) / ld.n_tnp;

    double p_a = (double)total_a / total_n;
    ld.p_b     = (double)ld.n_ena / total_n;
    ld.d       = ld_d(ena_a, ld.n_ena, total_a, total_n);

    // D'
    double d_max;
    if (ld.d > 0){
        d_max = std::min(p_a * (1 - ld.p_b), (1 - p_a) * ld.p_b);
    }
    else {
        d_max = std::min(p_a * ld.p_b, (1 - p_a) * (1 - ld.p_b));
    }
    if (d_max > 0){
        ld.d_prime  = ld.d / d_max;
        ld.dp_valid = true;
    }

    // r2
    double denom = p_a * (1 - p_a) * ld.p_b * (1 - ld.p_b);
    if (denom > 0){
        ld.r2       = ld.d * ld.d / denom;
        ld.r2_valid = true;
    }

    if (p_a == 0 || p_a == 1){
        return; // monomorphic in bear paw, D is 0 for every labeling
    }

    //
    // permutation test: keep the number of ENA mitotype fish (ena_k)
    // & relabel which fish they are. two sided on |D|
    //
    const double obs   = std::abs(ld.d);
    const double combs = n_choose_k(nfish, ena_k);
    uint extreme = 0;

    if (combs <= MAX_EXACT_PERMS){
        //
        // exact: every combination of ena_k fish out of nfish
        // idx holds the chosen fish, in increasing order
        //
        vector<uint> idx(ena_k);
        for (uint i = 0; i < (uint)ena_k; i++){
            idx[i] = i;
        }

        while (true){
            int group_a = 0, group_n = 0;
            for (uint i = 0; i < (uint)ena_k; i++){
                group_a += a_copies[idx[i]];
                group_n += n_copies[idx[i]];
            }
            double d = ld_d(group_a, group_n, total_a, total_n);
            if (std::abs(d) >= obs - PERM_TOL){
                extreme++;
            }
            ld.n_perm++;

            // move to the next combination
            int j = ena_k - 1;
            while (j >= 0 && idx[j] == nfish - ena_k + j){
                j--;
            }
            if (j < 0){
                break; // all combinations done
            }
            idx[j]++;
            for (uint i = j + 1; i < (uint)ena_k; i++){
                idx[i] = idx[i - 1] + 1;
            }
        }
        ld.p_perm = (double)extreme / ld.n_perm;
    }
    else {
        //
        // too many fish to enumerate, so shuffle the labels randomly
        // fixed seed so runs are reproducible
        //
        ld.exact = false;
        static std::mt19937 rng(12345);
        vector<uint> order(nfish);
        for (uint i = 0; i < nfish; i++){
            order[i] = i;
        }

        for (uint p = 0; p < N_RANDOM_PERMS; p++){
            std::shuffle(order.begin(), order.end(), rng);
            int group_a = 0, group_n = 0;
            for (uint i = 0; i < (uint)ena_k; i++){
                group_a += a_copies[order[i]];
                group_n += n_copies[order[i]];
            }
            double d = ld_d(group_a, group_n, total_a, total_n);
            if (std::abs(d) >= obs - PERM_TOL){
                extreme++;
            }
        }
        ld.n_perm = N_RANDOM_PERMS;
        ld.p_perm = (double)(extreme + 1) / (N_RANDOM_PERMS + 1);
    }
}

void
write_header(ofstream &fh, const vector<string> &pops){
    //
    // write the column headers for the output tsv
    //
    fh << "#Chrom\tBP\tGene_id\tRegion\tRef\tAlt\tClass\t"
       << pops[0] << "_allele\t"   << pops[1] << "_allele\t"
       << pops[0] << "_freq\t"     << pops[1] << "_freq\t"
       << pops[0] << "_nalleles\t" << pops[1] << "_nalleles";
    
    for (uint i = 0; i < 2; i++){
        string b = "BP_" + pops[i] + '_';
        fh << '\t'    << b << "homENA\t" << b
           << "het\t" << b << "homTNP\t" << b
           << "miss";
    }

    // for genotype summaries
    fh << '\t' << pops[0] << "_genotypes"
       << '\t' << "BP_"   << pops[0]      << "_genotypes"
       << '\t' << "BP_"   << pops[1]      << "_genotypes"
       << '\t' << pops[1] << "_genotypes" << '\n';

}

void
write_value(ofstream &fh, const bool valid, const double value){

    //
    // helper function to write a value, or NA if it is undefined
    //
    if (valid){
        fh << '\t' << value;
    }
    else {
        fh << "\tNA";
    }
}

void
write_ld_header(ofstream &fh, const vector<string> &pops){
    //
    // column headers for the mitonuclear LD output
    //
    fh << "#Chrom\tBP\tGene_id\tRegion\tRef\tAlt\tClass\tStatus\t"
       << pops[0] << "_allele\t" << pops[1]          << "_allele\t"
       << pops[0] << "_freq\t"   << pops[1]          << "_freq\t"
       << "BP_"   << pops[0]     << "_nalleles\tBP_" << pops[1] << "_nalleles\t"
       << "BP_"   << pops[0]     << "_q" << pops[0]  << "\tBP_" << pops[1] 
       << "_q"    << pops[0]     
       << "\tpB\tD\tD_prime\tr2\tn_perm\tperm_type\tp_perm\tBP_ndiploid\tBP_Ho\tBP_He\tBP_Fis";

    // for genotype summaries
    fh << '\t' << pops[0] << "_genotypes"
       << '\t' << "BP_"   << pops[0]      << "_genotypes"
       << '\t' << "BP_"   << pops[1]      << "_genotypes"
       << '\t' << pops[1] << "_genotypes" << '\n';

}

string
genotype_summary(const vector<Sample*> &samples, const vector<uint> &ploidy,
                 const vector<int> &alt_copies, const string &ref, const string &alt,
                 const mtType pop, const bool bear_paw){

    //
    // summary the genotypes of for the current pop (e.g., bear paw ENA)
    // entries will look something like "10 A/A, 2 A/T, 1 ./."
    //
    // counts:
    //  [0] ref/ref [1] ref/alt [2] alt/alt [3] ./.
    //  [4] ref     [5] alt     [6] .
    //

    uint counts[7] = {0, 0, 0, 0, 0, 0, 0};

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];

        if (sample == nullptr || ploidy[i] == 0){
            continue; // females at y do not contribute
        }

        if (sample->pop != pop || sample->bear_paw != bear_paw){
            continue; // not part of this group
        }

        if (ploidy[i] == 2){
            if (alt_copies[i] == -1){
                counts[3]++;
            }
            else {
                counts[alt_copies[i]]++;
            }
        }
        else {
            if (alt_copies[i] == -1){
                counts[6]++;
            }
            else {
                counts[4 + alt_copies[i]]++;
            }
        }
    } // end of i

    const string labels[7] = {ref + '/' + ref, ref + '/' + alt, alt + '/' + alt, "./.", ref, alt, "."};

    string summary;
    for (uint i = 0; i < 7; i++){
        if (counts[i] == 0){
            continue;
        }
        if (!summary.empty()){
            summary += ", ";
        }
        summary += to_string(counts[i]) + ' ' + labels[i];
    }

    if (summary.empty()){
        summary = "."; // no samples?
    }

    return summary;
}

void
write_site(ofstream &fh, const string &chrom, const uint pos, const Gene *gene,
           const RegionType region, const string &ref, const string &alt, const Classification c,
           const int ena_allele, const uint pure_alt[2], const uint pure_n[2],
           const uint bp[2][4], const vector<Sample*> &samples,
           const vector<uint> &ploidy, const vector<int> &alt_copies){
    //
    // write a passing site to the output tsv
    // freqs are of the ENA allele in pure ENA &
    // the TNP allele in pure TNP
    //

    double alt_freq_ena = (double)pure_alt[0] / pure_n[0];
    double alt_freq_tnp = (double)pure_alt[1] / pure_n[1];

    string ena_base, tnp_base;
    double ena_freq, tnp_freq;

    if (ena_allele == 1){
        ena_base = alt;
        tnp_base = ref;
        ena_freq = alt_freq_ena;
        tnp_freq = 1.0 - alt_freq_tnp;
    }
    else {
        ena_base = ref;
        tnp_base = alt;
        ena_freq = 1.0 - alt_freq_ena;
        tnp_freq = alt_freq_tnp;
    }

    if (c == Classification::BpHetSame){
        //
        // all non-bear paw samples are fixed
        // for the same allele
        //
        tnp_base = ena_base;
        if (ena_allele == 1){
            tnp_freq = alt_freq_tnp;
        }
        else {
            tnp_freq = 1.0 - alt_freq_tnp;
        }
    }

    string c_label;
    switch (c) {
        case Classification::Diff:
            c_label = "DIFF";
            break;
        case Classification::Private:
            c_label = "PRIVATE";
            break;
        case Classification::BpHetSame:
            c_label = "BP_HET_SAME";
            break;
        case Classification::BpHetDiff:
            c_label = "BP_HET_DIFF";
            break;
        default:
            c_label = "UNINFORMATIVE";
            break;
    }

    string gene_id = gene == nullptr ? "." : gene->id ;
  
    string r_label;
    switch (region) {
        case RegionType::Exonic:
            r_label = "exon";
            break;
        case RegionType::Promoter:
            r_label = "promoter";
            break;
        default:
            r_label = ".";
            break;
    }

    fh << chrom    << '\t' << pos      << '\t' << gene_id   << '\t'
       << r_label  << '\t' << ref      << '\t' << alt       << '\t'
       << c_label  << '\t'<< ena_base  << '\t' << tnp_base  << '\t'
       << ena_freq << '\t' << tnp_freq << '\t' << pure_n[0] << '\t' 
       << pure_n[1];

    for (uint i = 0; i < 2; i++){
        for (uint j = 0; j < 4; j++){
            fh << '\t' << bp[i][j];
        }
    }
    // genotype summaries: ENA, BP ENA, BP TNP, TNP
    fh << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::ENA, false)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::ENA, true)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::TNP, true)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::TNP, false)
       << '\n';

}

void
write_ld_site(ofstream &fh, const string &chrom, const uint pos, const Gene *gene,
              const RegionType region, const string &ref, const string &alt,
              const Classification c, const string &status, const int ena_allele,
              const uint pure_alt[2], const uint pure_n[2], const MitoLD &ld,
              const vector<Sample*> &samples, const vector<uint> &ploidy, const vector<int> &alt_copies){
    //
    // write one site to the mitonuclear LD output
    //

    string ena_base, tnp_base;
    double ena_freq, tnp_freq;

    if (ena_allele == 1){
        ena_base = alt;
        tnp_base = ref;
        ena_freq = (double)pure_alt[0] / pure_n[0];
        tnp_freq = 1.0 - (double)pure_alt[1] / pure_n[1];
    }
    else {
        ena_base = ref;
        tnp_base = alt;
        ena_freq = 1.0 - (double)pure_alt[0] / pure_n[0];
        tnp_freq = (double)pure_alt[1] / pure_n[1];
    }

    string c_label = c == Classification::Diff ? "DIFF" : "PRIVATE";
    string gene_id = gene == nullptr ? "." : gene->id;
    string r_label = ".";

    if (region == RegionType::Exonic){
        r_label = "exon";
    }
    else if (region == RegionType::Promoter){
        r_label = "promoter";
    }

    fh << chrom    << '\t' << pos      << '\t' << gene_id << '\t' << r_label << '\t'
       << ref      << '\t' << alt      << '\t' << c_label << '\t' << status  << '\t'
       << ena_base << '\t' << tnp_base << '\t' << ena_freq << '\t' << tnp_freq << '\t'
       << ld.n_ena << '\t' << ld.n_tnp;

    write_value(fh, ld.valid, ld.q_ena);
    write_value(fh, ld.valid, ld.q_tnp);
    write_value(fh, ld.valid, ld.p_b);
    write_value(fh, ld.valid, ld.d);
    write_value(fh, ld.dp_valid, ld.d_prime);
    write_value(fh, ld.r2_valid, ld.r2);

    if (ld.n_perm > 0){
        fh << '\t' << ld.n_perm;
        if (ld.exact){
            fh << "\texact";
        }
        else {
            fh << "\trandom";
        }
        fh << '\t' << ld.p_perm;
    }
    else {
        fh << "\t0\tNA\tNA";
    }

    fh << '\t' << ld.n_dip;
    write_value(fh, ld.n_dip > 0, ld.ho);
    write_value(fh, ld.n_dip > 0, ld.he);
    write_value(fh, ld.fis_valid, ld.fis);

    // genotype summaries: ENA, BP ENA, BP TNP, TNP
    fh << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::ENA, false)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::ENA, true)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::TNP, true)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::TNP, false)
       << '\n';
}


int
parse_vcf(const string &vcf, const double min_freq, const double min_frac, const double max_missing, const bool write_ld,
          const vector<string> &pops, unordered_map<string, Sample*> &pop_map, unordered_map<string, vector<Gene*>> &genome){

    //
    // this function will be the main work horse for this
    // code
    //
    //  1. allele frequencies in non-bear paw ENA & TNP samples
    //     define which allele belongs to ENA and which to B
    //  2. Bear Paw ENA samples must have a ENA allele (not homo TNP)
    //     & Bear Paw TNP samples must have a TNP allele (not homo ENA) 
    //

    bool gzipped = is_compressed(vcf);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, vcf);

    ofstream ofh("Bear_Paw_info_sites.tsv");

    if (!ofh.is_open()){
        cerr << "Error: Unable to create the output file in this directory\n";
        exit(EXIT_FAILURE);
    }

    ofh << std::fixed << setprecision(3);
    write_header(ofh, pops);

    // mitonuclear LD output if requested
    ofstream ld_fh;

    if (write_ld){
        ld_fh.open("Bear_Paw_mitonuclear_ld.tsv");
        if (!ld_fh.is_open()){
            cerr << "Error: Unable to create the LD output file in this directory\n";
            exit(EXIT_FAILURE);
        }
        ld_fh << std::fixed << setprecision(4);
        write_ld_header(ld_fh, pops);
    }

    MitoLD ld;
    uint   ld_sites = 0;
    
    //
    // sample -> pointer to Sample structs with pop & sex info
    //
    
    vector<string>  parts;
    vector<Gene*>  *genes;
    vector<Sample*> samples;
    vector<uint>    gene_ends, ploidy;
    vector<int>     alt_copies;
    Classification  clsfion;
    RegionType      region;
    string line, sample, chrom, curChrom;
    uint   pos, ngenes, furthest;
    bool   keep_reading;
    int    gt_idx, ena_allele;
    uint snp_sites = 0, informative = 0, shared_hom = 0, rescued = 0, het_diff = 0;
    uint passed = 0, miss_cnt = 0, exonic_sites = 0, promoter_sites = 0;
    uint pure_alt[2], pure_n[2]; // per site tallies for "pure" populations
    uint bp[2][4];               // per site bear paw genotype counts
    Gene *gene;

    while (true){
        keep_reading = read_line(gzipped, gz_fh, txt_fh, line);
        
        if (keep_reading == false){
            break;
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
                    if (pop_map.count(sample) == 0){
                        cerr << "Skipping " << sample << " in this analysis\n";
                        samples.push_back(nullptr);
                    }
                    else {
                        samples.push_back(pop_map[sample]);
                    }
                }
                ploidy.resize(samples.size());
                alt_copies.resize(samples.size());
            }
            continue;
        }

        if (samples.empty()){
            cerr << "Error: Could not find #CHROM header containg sample information\n";
            ofh.close();
            exit(1);
        }

        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }
        parse_tabular(line, parts);

        if (parts.size() != samples.size() + 9){
            cerr << "Error: Malformed vcf record. Expected " << samples.size() + 9
                 << " columns. Offending line:\n" << line << '\n';
            ofh.close();
            exit(EXIT_FAILURE);
        }

        chrom  = parts[0];
        pos    = (uint)stoi(parts[1]);
        gene   = nullptr;
        region = RegionType::NoRegion;

        // if using only genes, confirm that it is exonic about the exonic regions
        if (!genome.empty()){ 

            if (genome.count(chrom) == 0){
                continue; // no genes on this chrom
            }
            
            if (chrom != curChrom){
                // create the new vectors
                curChrom = chrom;
                genes    = &genome[chrom];
                ngenes   = genes->size();
                furthest = 0;
                gene_ends.clear();
                gene_ends.resize(ngenes);

                for (uint n = 0; n < ngenes; n++){
                    if ((*genes)[n]->search_end > furthest){
                        furthest = (*genes)[n]->search_end;
                    }
                    gene_ends[n] = furthest;
                }
            }
            
            // check if this position is exonic
            gene = is_candidate_region(pos, genes, gene_ends, region);
            if (gene == nullptr){
                continue;
            }
            if (region == RegionType::Exonic){
                exonic_sites++;
            }
            else if (region == RegionType::Promoter){
                promoter_sites++;
            }
        }
        
        // only focusing on biallelic sites
        const string &ref = parts[3];
        const string &alt = parts[4];

        if (!is_biallelic(ref, alt)){
            continue;
        }
        snp_sites++;

        gt_idx = find_gt_index(parts[8]);
        if (gt_idx == -1){
            continue; // no genotypes for this record
        }

        // call every sample & count the "pure" populations
        tally_genotypes(parts, samples, gt_idx, chrom, pos, ploidy, alt_copies, pure_alt, pure_n);

        // does this site have enough genotyped individuals
        if (!passes_missing(samples, ploidy, alt_copies, max_missing)){
            miss_cnt++;
            continue;
        }

        // which allele belongs to ENA. -1 if not informative
        clsfion    = Classification::Uinform; 
        ena_allele = classify_site(chrom, pure_alt, pure_n, min_freq, clsfion);

        //
        // mitonuclear LD for every informative nuclear site,
        // before the shared homozygote & bear paw filters.
        // the Status column says what those filters decided
        //
        if (write_ld && ena_allele != -1 && chrom != MT && chrom != Y){
            string status;
            if (shares_homozygotes(samples, ploidy, alt_copies)){
                status = "SHARED_HOM";
            }
            else {
                tally_bear_paw(samples, ploidy, alt_copies, ena_allele, bp);
                if (is_informative(bp, min_frac, clsfion, ena_allele, pure_alt, pure_n)){
                    status = "PASS";
                }
                else {
                    status = "BP_FAIL";
                }
            }

            mitonuclear_ld(samples, ploidy, alt_copies, ena_allele, ld);
            write_ld_site(ld_fh, chrom, pos, gene, region, ref, alt, clsfion, status, ena_allele,
                          pure_alt, pure_n, ld, samples, ploidy, alt_copies);
            ld_sites++;
        }

        // can we rescue this site?
        if (ena_allele != -1 && shares_homozygotes(samples, ploidy, alt_copies)){
            shared_hom++;
            ena_allele = -1; // no longer reportable as DIFF/PRIVATE
        }

        if (ena_allele == -1){
            //
            // we can attempt to rescue this hit if Bear Paw is
            // hetero while other pops are pure
            //
            if (chrom == MT || chrom == Y){
                continue; // haploid site
            }
            clsfion    = Classification::Uinform;
            ena_allele = classify_fixed_pure(samples, ploidy, alt_copies, clsfion);

            if (ena_allele == -1){
                continue; // a pure pop is not fixed
            }

            if (!bear_paw_all_het(samples, ploidy, alt_copies, min_frac)){
                continue;
            }

            rescued++;
            tally_bear_paw(samples, ploidy, alt_copies, ena_allele, bp);
            write_site(ofh, chrom, pos, gene, region, ref, alt, clsfion, ena_allele, pure_alt, pure_n, bp,
                       samples, ploidy, alt_copies);
            continue;
            
        }

        informative++;

        // does the bear paw samples carry their mtDNA's type?
        tally_bear_paw(samples, ploidy, alt_copies, ena_allele, bp);

        if (!is_informative(bp, min_frac, clsfion, ena_allele, pure_alt, pure_n)){
            continue;
        }
        passed++;

        //
        // a fixed difference where the bear paw samples are
        // heterozygous gets its own label (BP_HET_DIFF)
        //
        if (clsfion == Classification::Diff){
            Classification fixed_c = Classification::Uinform;
            int fixed_allele       = classify_fixed_pure(samples, ploidy, alt_copies, fixed_c);

            if (fixed_allele != -1 && fixed_c == Classification::BpHetDiff){
                if (bear_paw_all_het(samples, ploidy, alt_copies, min_frac)){
                    clsfion = Classification::BpHetDiff;
                    het_diff++;
                }
            }
        }

        write_site(ofh, chrom, pos, gene, region, ref, alt, clsfion, ena_allele, pure_alt, pure_n, bp, 
                   samples, ploidy, alt_copies);


    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);
    ofh.close();

    if (write_ld){
        ld_fh.close();
    }

    if (!genome.empty()){
        cerr << "Total exonic sites scanned "   << exonic_sites   << '\n'
             << "Total promoter sites scanned " << promoter_sites << '\n';
    }

    cerr << "Biallelic SNPs scanned "    << snp_sites << '\n'
         << "Sites removed for missing data " << miss_cnt << '\n'
         << "Sites removed for shared homozygotes " << shared_hom << '\n'
         << "Informative sites between " << pops[0]   << " and " << pops[1] << ' ' << informative << '\n'
         << "Sites passing the bear paw check " << passed << '\n'
         << "  of which fixed differences with heterozygous bear paw (BP_HET_DIFF) " << het_diff << '\n'
         << "Sites rescued (pure populations fixed, bear paw heterozygous) " << rescued << '\n';

    if (write_ld){
        cerr << "Sites written to the mitonuclear LD file " << ld_sites << '\n';
    }

    return 0;
}

void
help(){
    cerr << "Usage: ./find_bear_paw_hetero_sites -v vcf_file -p popmap.tsv -m <float> -f <float> -x <float> -u <int> -a ann.gtf [optional] -g gene_list.txt [optional] -l [optional]\n"
         << "\t-v: vcf file containing the bear paw samples along with the other TNP and ENA stickleback population [required]\n"
         << "\t-p: A population map containing 3 columns: [Sample][Population][Sex] [required]\n"
         << "\t-m: Minimum frequency of an allele in a population (0.5 < x <= 1), to consider it differentiated [default: 1.0]\n"
         << "\t-f: Minimum fraction of Bear Paw Samples to carry the specific mt-Type allele [default: 1.0]\n"
         << "\t-x: Maximum fraction of missing samples allowed in each group [default: 1.0; no filter]\n"
         << "\t-a: The genome annotation file to overlap SNPs with\n"
         << "\t-g: A file contain the gene_ids (one per line) of interest. Only these genes will be investigated in the annotation\n"
         << "\t-u: Number of sites to look upstream of a gene when using a genome annotation [default: 0]\n"
         << "\t-l: Also write Bear_Paw_mitonuclear_ld.tsv: LD between the nuclear ENA allele and the ENA mitotype within Bear Paw for every informative (DIFF/PRIVATE) nuclear site, with a permutation p-value [default: off]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, pop_file, arg, genes_list;
    double min_freq = 1.0, min_frac = 1.0, max_missing = 1.0;
    int  upstream = 0;
    bool write_ld = false;

    // expect at least a single argument
    if (argc < 2){
        help();
    }

    for (int i = 1; i < argc; i++){
        arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = argv[i + 1];
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = argv[i + 1];
        }
        else if (arg == "-p" && i + 1 < argc){
            pop_file = argv[i + 1];
        }
        else if (arg == "-g" && i + 1 < argc){
            genes_list = argv[i + 1];
        }
        else if (arg == "-l"){
            write_ld = true;
        }
        else if (arg == "-h"){
            help();
        }
        else if (arg == "-m" && i + 1 < argc){
            try {
                min_freq = stod(argv[i + 1]);

                if (min_freq <= 0.5 || min_freq > 1.0){
                    cerr << "Error: -m must be a value greater than 0.5 and at most 1.0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::invalid_argument &e){
                cerr << "Error: -m must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
        }
        else if (arg == "-f" && i + 1 < argc){
            try {
                min_frac = stod(argv[i + 1]);

                if (min_frac < 0 || min_frac > 1.0){
                    cerr << "Error: -f must be a value within 0 and 1.0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::invalid_argument &e){
                cerr << "Error: -f must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
        }
        else if (arg == "-x" && i + 1 < argc){
            try {
                max_missing = stod(argv[i + 1]);

                if (max_missing < 0 || max_missing > 1.0){
                    cerr << "Error: -x must be a value within 0 and 1.0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::invalid_argument &e){
                cerr << "Error: -x must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
        }
        else if (arg == "-u" && i + 1 < argc){
            try {
                upstream = stoi(argv[i + 1]);

                if (upstream < 0){
                    cerr << "Error: -u cannot be less than 0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::invalid_argument &e){
                cerr << "Error: -u must be an integer. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
        }
    }
    if (vcf.empty()){
        cerr << "-v vcf_file is required\n";
        exit(EXIT_FAILURE);
    }
    else if (!file_exists(vcf)){
        cerr << "Unable to find " << vcf << '\n';
        exit(1);
    }
    if (pop_file.empty()){
        cerr << "-p pop_file is required\n";
        exit(EXIT_FAILURE);
    }
    else if (!file_exists(pop_file)){
        cerr << "Unable to find " << pop_file << '\n';
        exit(1);
    }

    //
    // check if we will use an annotation file
    // to only look at exonic regions for
    // protein coding genes
    //
    unordered_map<string, vector<Gene*>> genome;
    bool use_ann = !ann.empty();
    if (use_ann){
        if (!file_exists(ann)){
            cerr << "Unable to find " << ann << '\n';
            exit(1);
        }

        bool only_targets = !genes_list.empty();
        // only use a subset of genes based on their gene_ids
        if (only_targets){
            if (!file_exists(genes_list)){
                cerr << "Unable to find " << genes_list << '\n';
                exit(EXIT_FAILURE);
            }
        }
        unordered_set<string> target_genes;
        load_target_genes(genes_list, only_targets, target_genes);
        parse_annotation(ann, genome, only_targets, target_genes, (uint)upstream);
    }
    else if (upstream > 0){
        cerr << "Ignoring -u parameter\n";
    }

    // load the samples
    unordered_map<string, Sample*> pop_map;
    vector<string> pops; // index 0 -> ENA, index 1 -> TNP
    parse_popmap(pop_file, pop_map, pops);
    
    // parse the vcf file
    parse_vcf(vcf, min_freq, min_frac, max_missing, write_ld, pops, pop_map, genome);

    // we are done
    clean_up(genome, pop_map);

    return 0;
}