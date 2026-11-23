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
#include <utility>
#include <cmath>
#include <random>
#include <map>

using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::map;
using std::pair;
using std::cerr;
using std::cout;
using std::stoi;
using std::sort;
using std::find;
using std::make_pair;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::setprecision;
using std::stod;
using std::stoi;

//
// this code is specifically designed. for finding heterozygous
// INDELS for the population of interest (bear paw). If given an
// annotation file, only exonic sides will be searched. If given a
// list of genes, then only those genes within the annotation will be
// searched
//

//
// records are collapsed if they overlap. Alleles are compared by
// the sequence they create across the whole locus, so the
// same event written by two records is the same allele
//

typedef unsigned int uint;

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
    BpHetDiff = 1, // pure pops fixed for different alleles, bear paw het
    BpHet     = 2, // recused
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

//
// structs for handling loci spanning >= 1 record
// Record will be for a single vcf record
// Locus will handle loci spanning > 1 record
//
struct Record {
    uint pos;
    uint end;
    string ref;
    vector<string> alts;
    vector<string> parts;
};

struct Locus {
    string         chrom;
    uint           start;
    uint           end;
    uint           nrecords;
    vector<string> alleles; // [0] = ref, [k] = alt seq
    vector<uint>   ploidy;  // ordered by sample
    vector<int>    a1;      // alleles follow sample order, -1 = missing
    vector<int>    a2; 
    vector<double> ab;      // allele balance, -1 not hetero or no depth
};

//
// options from the user to keep
// things organized 
//
struct Params {
    double min_freq    = 1.0; // -m
    double min_frac    = 1.0; // -f
    double max_missing = 1.0; // -x
    int    min_gq      = 0;   // -q, 0 = off
    double min_ab      = 0;   // -b, 0 = off
};

//
// keep a singleton struct that holds
// all the counters to report at the end
//
struct Counter {
    uint records      = 0;
    uint loci         = 0;
    uint multi_record = 0; // loci built from more than one record
    uint haploid_skip = 0; // MT & Y loci
    uint outside      = 0; // loci outside the annotation regions
    uint exonic       = 0;
    uint promoter     = 0;
    uint ref_mismatch = 0;
    uint no_alleles   = 0; // only spanning deletions (*) in the locus
    uint conflicts    = 0; // sample calls that could not be combined
    uint low_gq       = 0; // sample calls removed by -q
    uint low_ab       = 0; // sample het calls removed by -b
    uint too_missing  = 0;
    uint diff         = 0; // pure populations differ (-m)
    uint shared_hom   = 0; // of those, both pure pops share a homozygote
    uint bp_het_diff  = 0;
    uint bp_het       = 0;
};

//
// for Bear Paw samples only
// index 0 = ENA; index 1 = TNP
//
struct BpStats {
    uint   called[2]   = {0, 0}; // called samples (any ploidy)
    uint   diploid[2]  = {0, 0}; // called diploid samples
    uint   het[2]      = {0, 0}; // het for any two alleles
    uint   het_diff[2] = {0, 0}; // het for exactly the ENA & TNP alleles
    uint   miss[2]     = {0, 0};
    uint   ab_n        = 0;      // het calls with allele depths
    uint   ab_low      = 0;      // of those, allele balance below LOW_AB
    double ab_sum      = 0;
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

    bool is_exonic(const uint start, const uint end){
        //
        // exons are ordered,
        // so just see if there is an overlap
        //

        for (uint i = 0; i < this->exons.size(); i++){
            const Exon &e = this->exons[i];
            if (e.start <= end && start <= e.end){
                return true;
            }
        }
        // no overlap
        return false;
    }

    bool in_promoter(const uint start, const uint end){
        //
        // is this site in the promoter search space
        //
        if (this->has_promoter){
            return this->prom_start <= end && start <= this->prom_end;
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
const double LOW_AB    = 0.2;        // het calls below this allele balance are reported as low

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
parse_string(const string &s, const char sep, vector<string> &out){

    //
    // split s on sep. an empty string gives an empty vector
    //

    out.clear();
    if (s.empty()){
        return 1;
    }

    uint start = 0;
    for (uint i = 0; i <= s.size(); i++){
        if (i == s.size() || s[i] == sep){
            out.push_back(s.substr(start, i - start));
            start = i + 1;
        }
    }

    return 0;
}

int
to_int(const string &s){

    //
    // convert a vcf number to an int
    // missing ('.' or empty) or invalid returns -1
    //

    if (s.empty() || s == "."){
        return -1;
    }
    for (uint i = 0; i < s.size(); i++){
        // compare ascii values
        if (s[i] < '0' || s[i] > '9'){
            return -1;
        }
    }
    return stoi(s);
}

bool
parse_gt(const string &gt, vector<int> &alleles){

    //
    // collect the allele indices of a GT subfield (e.g. 0/1, 1|0, 1)
    // return false if any is missing ('.') or invalid
    //

    alleles.clear();
    if (gt.empty()){
        return false;
    }

    string allele;
    for (uint i = 0; i <= gt.size(); i++){
        if (i == gt.size() || gt[i] == '/' || gt[i] == '|'){
            int a = to_int(allele);
            if (a < 0){
                return false; // missing or malformed
            }
            alleles.push_back(a);
            allele.clear();
        }
        else {
            allele += gt[i];
        }
    }

    return true;
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

int
set_chrom_genes(const string &chrom, unordered_map<string, vector<Gene*>> &genome,
                vector<Gene*> *&genes, vector<uint> &gene_ends){

    //
    // helper func() to populate the genes_end array
    // in order to use binary search to find genes extending
    // long distances
    //

    genes = &genome[chrom];

    const uint ngenes = genes->size();
    uint furthest     = 0;

    gene_ends.clear();
    gene_ends.resize(ngenes);

    for (uint n = 0; n < ngenes; n++){
        if ((*genes)[n]->search_end > furthest){
            furthest = (*genes)[n]->search_end;
        }
        gene_ends[n] = furthest;
    }

    return 0;
}

Gene *
is_candidate_region(const uint start, const uint end, vector<Gene*>*genes, vector<uint> &gene_ends, RegionType &region){

    //
    // helper function to determine
    // whether this locus overlaps exonic region, a candidate promoter region, 
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
        if (gene_ends[mid] >= start){
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
        
        if (end < gene->search_start){
            break; // this & all later genes begin after the locus
        }
        if (start <= gene->search_end){
            if (gene->start <= end && start <= gene->end && gene->is_exonic(start, end)){
                region = RegionType::Exonic;
                return gene; // lands in an exon
            }
            else if (promoter_hit == nullptr && gene->in_promoter(start, end)){
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

bool
passes_missing(const vector<Sample*> &samples, const Locus &locus, const double max_missing){
    //
    // check that each group has at most max_missing
    // of its samples missing at this locus
    //
    // expected[pop][bp] & missing[pop][bp]
    // pop: 0 = ENA, 1 = TNP. bp: 0 = pure, 1 = bear paw
    //

    uint expected[2][2] = {{0, 0}, {0, 0}};
    uint missing[2][2]  = {{0, 0}, {0, 0}};

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || locus.ploidy[i] == 0){
            continue;
        }

        uint bp = 0;
        if (sample->bear_paw){
            bp = 1;
        }

        expected[sample->pop][bp]++;
        if (locus.a1[i] < 0){
            missing[sample->pop][bp]++;
        }
    }

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

    return true;
}

int
find_format_index(const string &format, const string &key){

    //
    // find which ':' delimited subfield
    // holds the key (e.g., GT, AD, DP, GQ)
    //

    uint start = 0, idx = 0, length = format.size();

    for (uint i = 0; i <= length; i++){
        if (i == length || format[i] == ':'){
            int len  = i - start; 
            if (len == key.size() && format.compare(start, i - start, key) == 0){
                return idx;
            }
            idx++;
            start = i + 1;
        }
    }

    // key not found
    return -1;
}

string
get_subfield(const string &field, const int idx){

    //
    // return the string at the ':' delimited subfield of a sample
    // column, or an empty string if it is not there
    //

    if (idx < 0){
        return "";
    }

    int start = 0, n = 0, length = field.size();

    for (int i = 0; i <= length; i++){
        if (i == length || field[i] == ':'){
            if (n == idx){
                return field.substr(start, i - start);
            }
            n++;
            start = i + 1;
        }
    }

    return "";
}

bool
build_locus(const vector<Record> &records, const string &chrom, const vector<Sample*> &samples,
            const Params &params, Locus &locus, Counter &counts){

    //
    // collapse the overlapping records into one locus
    //
    //  1. rebuild the reference sequence across the locus from the REF
    //     fields (they overlap, so they cover the locus without gaps)
    //  2. turn every alt into the sequence it creates across the locus,
    //     identical sequences are the same allele
    //  3. combine each sample's calls across the records into one genotype
    //
    // returns false if the locus can't be built
    //

    locus.chrom    = chrom;
    locus.start    = records[0].pos;
    locus.end      = records[0].end;
    locus.nrecords = records.size();

    // this will be used repetively
    const uint nrec = records.size();

    for (uint r = 1; r < nrec; r++){
        locus.end = locus.end > records[r].end ? locus.end : records[r].end;
    }

    //
    // 1. reference sequence across the locus
    //
    const uint length = locus.end - locus.start + 1;
    string ref_seq(length, 'N');
    vector<bool> filled(length, false);

    for (uint r = 0; r < nrec; r++){
        const Record &rec = records[r];
        uint offset       = rec.pos - locus.start;

        for (uint j = 0; j < rec.ref.size(); j++){
            char base = toupper(rec.ref[j]);

            //
            // the nucleotide was been filled, so we can't be 100% sure
            // even a single bp is considered an indel
            //
            if (filled[offset + j] && ref_seq[offset + j] != base){
                cerr << "Warning: REF bases disagree at " << chrom << ':' << rec.pos
                     << ", skipping this locus\n";
                counts.ref_mismatch++;
                return false;
            }
            
            ref_seq[offset + j] = base;
            filled[offset + j]  = true;
        }
    }

    //
    // 2. alleles as sequences across the locus
    // allele_idx[r][g] = locus allele for GT index g of record r
    // -2 = spanning deletion (*), the copy is described by another record
    //
    locus.alleles.clear();
    locus.alleles.push_back(ref_seq);

    unordered_map<string,int> seq_to_idx;
    seq_to_idx[ref_seq] = 0; // map of allele to index

    vector<vector<int>> allele_idx(nrec);

    for (uint r = 0; r < nrec; r++){
        const Record &rec = records[r];
        uint offset       = rec.pos - locus.start;

        allele_idx[r].push_back(0); // GT index 0 is the reference

        for (uint j = 0; j < rec.alts.size(); j++){
            const string &alt = rec.alts[j];
            if (alt == "*" || alt == "."){ // missing
                allele_idx[r].push_back(-2);
                continue;
            }

            string upper_alt = alt;
            for (uint k = 0; k < upper_alt.size(); k++){
                upper_alt[k] = toupper(upper_alt[k]);
            }

            string seq = ref_seq.substr(0, offset) + upper_alt + ref_seq.substr(offset + rec.ref.size());

            if (seq_to_idx.count(seq) == 0){
                seq_to_idx[seq] = locus.alleles.size();
                locus.alleles.push_back(seq);
            }
            allele_idx[r].push_back(seq_to_idx[seq]);
        }
    }

    // ref will occupy 1 space
    if (locus.alleles.size() < 2){
        counts.no_alleles++;
        return false; // nothing but spanning deletions
    }

    //
    // where are GT, AD, DP & GQ in each record
    //
    vector<int> gt_idx(nrec), ad_idx(nrec), dp_idx(nrec), gq_idx(nrec);

    for (uint r = 0; r < nrec; r++){
        const string &format = records[r].parts[8];
        gt_idx[r] = find_format_index(format, "GT");
        ad_idx[r] = find_format_index(format, "AD");
        dp_idx[r] = find_format_index(format, "DP");
        gq_idx[r] = find_format_index(format, "GQ");
    }

    //
    // 3. one genotype per sample
    //
    const uint nsamples = samples.size();
    const uint nalleles = locus.alleles.size();

    locus.ploidy.assign(nsamples, 0);
    locus.a1.assign(nsamples, -1);
    locus.a2.assign(nsamples, -1);
    locus.ab.assign(nsamples, -1);

    vector<int>    gt;
    vector<string> ad_parts;
    vector<int>    copies(nalleles), rec_copies(nalleles), depth(nalleles);

    for (uint i = 0; i < nsamples; i++){
        Sample *sample = samples[i];
        if (sample == nullptr){
            continue; // not in pop map
        }

        const int ploidy = expected_ploidy(chrom, locus.start, sample->sex);
        locus.ploidy[i]  = ploidy;
        if (ploidy == 0){
            continue; // no call expected (females on the Y)
        }

        copies.assign(nalleles, 0);
        depth.assign(nalleles, 0);
        bool missing    = false;
        bool ref_depth  = false; // has the reference depth been set
        int  first_ad   = -1;    // reference depth of the first record with depths

        for (uint r = 0; r < nrec && !missing; r++){
            const string &field = records[r].parts[i + 9];

            if (!parse_gt(get_subfield(field, gt_idx[r]), gt)){
                missing = true;
                break;
            }

            // a call with no reads is not a call
            int dp = to_int(get_subfield(field, dp_idx[r]));
            if (dp == 0){
                missing = true;
                break;
            }

            // does this pass the genotype quality scoring
            if (params.min_gq > 0){
                int gq = to_int(get_subfield(field, gq_idx[r]));
                if (gq >= 0 && gq < params.min_gq){
                    counts.low_gq++;
                    missing = true;
                    break;
                }
            }

            // haploid regions may be written as 1 or 1/1
            if (ploidy == 1 && gt.size() == 2){
                if (gt[0] != gt[1]){
                    missing = true; // het call in a haploid region
                    break;
                }
                gt.pop_back();
            }
            if ((int)gt.size() != ploidy){
                missing = true; // unexpected ploidy
                break;
            }

            //
            // copies of each locus allele in this record. a record
            // that repeats an allele from another record (same sequence)
            // should not count it twice, so keep the max across records
            //
            rec_copies.assign(nalleles, 0);
            bool carries_alt = false;

            for (uint j = 0; j < gt.size(); j++){
                if (gt[j] >= (int)allele_idx[r].size()){
                    missing = true; // GT index without an alt
                    break;
                }
                int k = allele_idx[r][gt[j]];
                if (k <= 0){
                    continue; // reference or spanning deletion
                }
                rec_copies[k]++;
                carries_alt = true;
            }
            if (missing){
                break;
            }

            for (uint k = 1; k < nalleles; k++){
                if (rec_copies[k] > copies[k]){
                    copies[k] = rec_copies[k];
                }
            }

            //
            // allele depths: AD[0] is the reference, AD[g] the g'th alt
            // the reference depth comes from the record where the sample
            // carries an alt, or the first record if it never does
            //
            parse_string(get_subfield(field, ad_idx[r]), ',', ad_parts);

            if (ad_parts.size() == allele_idx[r].size()){
                int ref_ad = to_int(ad_parts[0]);

                // likely missing, so default to ref
                if (first_ad == -1){
                    first_ad = ref_ad;
                
                }
                if (carries_alt && !ref_depth && ref_ad >= 0){
                    depth[0]  = ref_ad;
                    ref_depth = true;
                }

                for (uint g = 1; g < ad_parts.size(); g++){
                    int k  = allele_idx[r][g];
                    int ad = to_int(ad_parts[g]);
                    if (k > 0 && ad > depth[k]){
                        depth[k] = ad;
                    }
                }
            }
        } // end of r

        if (missing){
            continue;
        }

        if (!ref_depth && first_ad >= 0){
            depth[0] = first_ad;
        }

        // how many non reference copies across the locus
        int total = 0;
        for (uint k = 1; k < nalleles; k++){
            total += copies[k];
        }

        if (total > ploidy){
            counts.conflicts++; // e.g. 0/1 in three overlapping records
            continue;
        }

        // build the genotype, reference copies fill the rest
        vector<int> geno;
        for (uint k = 1; k < nalleles; k++){
            for (int c = 0; c < copies[k]; c++){
                geno.push_back(k);
            }
        }
        while ((int)geno.size() < ploidy){
            geno.push_back(0);
        }
        sort(geno.begin(), geno.end());

        if (ploidy == 1){
            locus.a1[i] = geno[0];
            continue;
        }

        // het calls get an allele balance (minor allele fraction)
        double ab = -1;
        if (geno[0] != geno[1]){
            int d1 = depth[geno[0]];
            int d2 = depth[geno[1]];
            if (d1 + d2 > 0){
                if (d1 < d2){
                    ab = (double)d1 / (d1 + d2);
                }
                else {
                    ab = (double)d2 / (d1 + d2);
                }
            }
            if (params.min_ab > 0 && ab >= 0 && ab < params.min_ab){
                counts.low_ab++;
                continue; // treated as missing
            }
        }

        locus.a1[i] = geno[0];
        locus.a2[i] = geno[1];
        locus.ab[i] = ab;
    } // end of i

    return true;
}

int
tally_pure(const vector<Sample*> &samples, const Locus &locus,
           vector<uint> pure_counts[2], uint pure_n[2]){

    //
    // count the copies of every allele in the pure populations
    // pure_counts[pop][allele] & pure_n[pop] (allele copies called)
    // index 0 = ENA, 1 = TNP
    //

    const uint nalleles = locus.alleles.size();

    for (uint p = 0; p < 2; p++){
        pure_counts[p].assign(nalleles, 0);
        pure_n[p] = 0;
    }

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || sample->bear_paw || locus.a1[i] < 0){
            continue; // only called pure samples
        }

        pure_counts[sample->pop][locus.a1[i]]++;
        pure_n[sample->pop]++;

        if (locus.a2[i] >= 0){
            pure_counts[sample->pop][locus.a2[i]]++;
            pure_n[sample->pop]++;
        }
    }

    return 0;
}

int
tally_bear_paw(const vector<Sample*> &samples, const Locus &locus,
               const int ena_allele, const int tnp_allele, BpStats &bp){

    //
    // count the bear paw genotypes in each group (BP ENA, BP TNP)
    // het_diff counts hets for exactly the ENA & TNP alleles
    // only diploid calls can be heterozygous
    //

    bp = BpStats(); // reset
    uint idx;

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || !sample->bear_paw || locus.ploidy[i] == 0){
            continue;
        }

        idx = sample->pop;
        if (locus.a1[i] < 0){
            bp.miss[idx]++;
            continue;
        }
        bp.called[idx]++;

        if (locus.a2[i] < 0){
            continue; // haploid
        }
        bp.diploid[idx]++;

        if (locus.a1[i] == locus.a2[i]){
            continue; // homozygous
        }
        bp.het[idx]++;

        bool ena_tnp = (locus.a1[i] == ena_allele && locus.a2[i] == tnp_allele) ||
                       (locus.a1[i] == tnp_allele && locus.a2[i] == ena_allele);
        if (ena_tnp){
            bp.het_diff[idx]++;
        }

        if (locus.ab[i] >= 0){
            bp.ab_n++;
            bp.ab_sum += locus.ab[i];
            if (locus.ab[i] < LOW_AB){
                bp.ab_low++;
            }
        }
    }

    return 0;
}

int
major_allele(const vector<uint> &counts){

    //
    // the most common allele (lowest index wins a tie)
    // return -1 if nothing was called
    //

    int  best   = -1;
    uint best_n = 0;

    for (uint i = 0; i < counts.size(); i++){
        if (counts[i] > best_n){
            best   = i;
            best_n = counts[i];
        }
    }

    return best;
}

bool
classify_diff(const vector<uint> pure_counts[2], const uint pure_n[2], const double min_freq,
              int &ena_allele, int &tnp_allele){

    //
    // determine if the the pure populations differ at this locus
    // the ENA allele must be at >= min_freq in pure ENA & <= 1 - min_freq
    // in pure TNP, & the TNP allele the reverse. min_freq > 0.5, so there
    // is at most one of each
    //

    ena_allele = major_allele(pure_counts[0]);
    tnp_allele = major_allele(pure_counts[1]);

    if (pure_n[0] == 0 || pure_n[1] == 0){
        return false; // no calls in one of the "pure" populations
    }

    int ena = -1, tnp = -1, nalleles = pure_counts[0].size();

    for (int i = 0; i < nalleles; i++){
        double freq_ena = (double)pure_counts[0][i] / pure_n[0];
        double freq_tnp = (double)pure_counts[1][i] / pure_n[1];

        if (freq_ena >= min_freq && freq_tnp <= 1.0 - min_freq){
            ena = i;
        }
        if (freq_tnp >= min_freq && freq_ena <= 1.0 - min_freq){
            tnp = i;
        }
    }

    // neither was set
    if (ena == -1 || tnp == -1){
        return false;
    }

    ena_allele = ena;
    tnp_allele = tnp;
    return true;
}

bool
bear_paw_het_diff(const BpStats &bp, const double min_frac){
    //
    // at least min_frac of the called diploid bear paw samples
    // in each group must be het for the ENA & TNP alleles
    //

    for (uint g = 0; g < 2; g++){
        if (bp.diploid[g] == 0){
            return false; // can't evaluate
        }
        double frac = (double)bp.het_diff[g] / bp.diploid[g];
        if (frac < min_frac){
            return false;
        }
    }
    return true;
}

bool
shares_homozygotes(const vector<Sample*> &samples, const Locus &locus){

    //
    // are there homozygous individuals for the same allele in both pure populations?
    //
    // haploid calls (male X) count as homozygous
    //
    // hom is really a 2D vector, with the following mapping
    // hom[pop][allele], pop: 0 = ENA, 1 = TNP
    //

    const uint nalleles = locus.alleles.size();
    vector<bool> hom[2];
    hom[0].assign(nalleles, false); // resize each vec by the num of alleles
    hom[1].assign(nalleles, false);

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || sample->bear_paw || locus.a1[i] < 0){
            continue; // only called pure samples
        }
        // alleles sorted by samples
        if (locus.a2[i] < 0 || locus.a1[i] == locus.a2[i]){
            hom[sample->pop][locus.a1[i]] = true;
        }
    }

    for (uint k = 0; k < nalleles; k++){
        if (hom[0][k] && hom[1][k]){
            return true;
        }
    }
    return false;
}

int
write_header(ofstream &fh, const vector<string> &pops){
    //
    // write the column headers for the output tsv
    //
    fh << "#Chrom\tStart\tEnd\tGene_id\tRegion\tN_records\tAlleles\tClass\tPure_class\t"
       << pops[0] << "_allele\t"   << pops[1] << "_allele\t"
       << pops[0] << "_freq\t"     << pops[1] << "_freq\t"
       << pops[0] << "_nalleles\t" << pops[1] << "_nalleles";

    for (uint i = 0; i < 2; i++){
        string b = "BP_" + pops[i] + '_';
        fh << '\t' << b     << "called\t" << b       << "het\t"
           << b    << "het" << pops[0]    << pops[1] << '\t' 
           << b    << "miss";
    }

    fh << "\tBP_het_meanAB\tBP_het_lowAB";

    // genotype summaries for each group
    fh << '\t' << pops[0] << "_genotypes"
       << '\t' << "BP_"   << pops[0]      << "_genotypes"
       << '\t' << "BP_"   << pops[1]      << "_genotypes"
       << '\t' << pops[1] << "_genotypes" << "\tAllele_seqs\n";

    return 0;
}

string
allele_summary(const Locus &locus, const RegionType region){

    //
    // describe each alt by its length change relative to the reference (e.g. "1:-22,2:+34")
    // Exonic loci also get fs (frameshift, length change not a multiple of 3) or if (in frame)
    //

    string    summary;
    const int ref_len = locus.alleles[0].size();

    for (uint k = 1; k < locus.alleles.size(); k++){
        int change = (int)locus.alleles[k].size() - ref_len;

        if (!summary.empty()){
            summary += ',';
        }
        summary += to_string(k) + ':';
        if (change > 0){
            summary += '+';
        }
        summary += to_string(change);

        if (region == RegionType::Exonic){
            if (change % 3 != 0){
                summary += ":fs"; // frameshift
            }
            else {
                summary += ":if"; // in frame
            }
        }
    }

    return summary;
}

string
genotype_summary(const vector<Sample*> &samples, const Locus &locus,
                 const mtType pop, const bool bear_paw){
    //
    // summarize the genotypes of one group 
    // Alleles column, e.g. "10 0/0, 2 0/1, 1 1/2, 1 ./."
    //

    // use a map b/c we care about orientation order
    map<pair<int,int>, uint> counts; // (a1, a2) -> count, a2 = -1 for haploid
    uint miss_dip = 0, miss_hap = 0;

    for (uint i = 0; i < samples.size(); i++){
        Sample *sample = samples[i];
        if (sample == nullptr || locus.ploidy[i] == 0){
            continue;
        }
        else if (sample->pop != pop || sample->bear_paw != bear_paw){
            continue; // not in this group
        }

        if (locus.a1[i] < 0){
            if (locus.ploidy[i] == 2){
                miss_dip++;
            }
            else {
                miss_hap++;
            }
            continue;
        }
        // genotype num, count
        pair<int,int> p = make_pair(locus.a1[i], locus.a2[i]);
        counts[p]++;
    }

    string summary;
    for (auto itr = counts.begin(); itr != counts.end(); itr++){
        if (!summary.empty()){
            summary += ", ";
        }
        summary += to_string(itr->second) + ' ' + to_string(itr->first.first);
        if (itr->first.second >= 0){
            summary += '/' + to_string(itr->first.second);
        }
    }

    if (miss_dip > 0){
        if (!summary.empty()){
            summary += ", ";
        }
        summary += to_string(miss_dip) + " ./.";
    }
    if (miss_hap > 0){
        if (!summary.empty()){
            summary += ", ";
        }
        summary += to_string(miss_hap) + " .";
    }

    if (summary.empty()){
        summary = "."; // no samples expected in this group
    }

    return summary;
}

void
write_locus(ofstream &fh, const Locus &locus, const Gene *gene, const RegionType region,
            const Classification c, const string &pure_class, const int ena_allele,
            const int tnp_allele, const vector<uint> pure_counts[2], const uint pure_n[2],
            const BpStats &bp, const vector<Sample*> &samples){
    //
    // write a locus to the output tsv
    // freqs are of the ENA allele in pure ENA &
    // the TNP allele in pure TNP
    //

    string c_label, r_label;
    switch (c){
        case Classification::BpHetDiff:
            c_label = "BP_HET_DIFF";
            break;
        case Classification::BpHet:
            c_label = "BP_HET";
            break;
        default:
            c_label = "UNINFORMATIVE"; // should not happen
            break;
    }

    switch (region){
        case RegionType::Exonic:
            r_label = "exon";
            break;
        case RegionType::Promoter:
            r_label = "promoter";
            break;
        default:
            r_label = "."; // no annotated used
            break;
    }

    string gene_id = gene == nullptr ?  "." : gene->id;


    fh << locus.chrom << '\t' << locus.start << '\t' << locus.end      << '\t'
       << gene_id     << '\t' << r_label     << '\t' << locus.nrecords << '\t'
       << allele_summary(locus, region)      << '\t' << c_label        << '\t' 
       << pure_class;

    // the major allele of each pure population (or "." if none called)
    if (ena_allele >= 0){
        fh << '\t' << ena_allele;
    }
    else {
        fh << "\t.";
    }

    if (tnp_allele >= 0){
        fh << '\t' << tnp_allele;
    }
    else {
        fh << "\t.";
    }

    // ena frequency
    if (ena_allele >= 0 && pure_n[0] > 0){
        fh << '\t' << (double)pure_counts[0][ena_allele] / pure_n[0];
    }
    else {
        fh << "\tNA";
    }

    // tnp frequency
    if (tnp_allele >= 0 && pure_n[1] > 0){
        fh << '\t' << (double)pure_counts[1][tnp_allele] / pure_n[1];
    }
    else {
        fh << "\tNA";
    }

    // total counts
    fh << '\t' << pure_n[0] << '\t' << pure_n[1];

    for (uint g = 0; g < 2; g++){
        fh << '\t' << bp.called[g] << '\t' << bp.het[g] << '\t'
           << bp.het_diff[g] << '\t' << bp.miss[g];
    }

    if (bp.ab_n > 0){
        fh << '\t' << bp.ab_sum / bp.ab_n;
    }
    else {
        fh << "\tNA";
    }
    fh << '\t' << bp.ab_low;

    // genotype summaries: ENA, BP ENA, BP TNP, TNP
    fh << '\t' << genotype_summary(samples, locus, mtType::ENA, false)
       << '\t' << genotype_summary(samples, locus, mtType::ENA, true)
       << '\t' << genotype_summary(samples, locus, mtType::TNP, true)
       << '\t' << genotype_summary(samples, locus, mtType::TNP, false)
       << '\t';

    // the allele sequences across the locus, reference first
    for (uint k = 0; k < locus.alleles.size(); k++){
        if (k > 0){
            fh << ',';
        }
        fh << locus.alleles[k];
    }

    fh << '\n';
}

int
process_records(const vector<Record> &records, const string &chrom,
                unordered_map<string, vector<Gene*>> &genome, vector<Gene*> *&genes,
                vector<uint> &gene_ends, string &gene_chrom, const vector<Sample*> &samples,
                const Params &params, ofstream &ofh, Counter &counts){

    //
    // collapse a collections of overlapping records & report it
    // if the bear paw samples are heterozygous
    //

    if (records.empty()){
        return 1;
    }

    counts.loci++;
    if (records.size() > 1){
        counts.multi_record++;
    }

    if (chrom == MT || chrom == Y){
        counts.haploid_skip++;
        return 1; // haploid, no heterozygotes possible
    }

    // locus span, needed for the annotation check
    uint start = records[0].pos, end = records[0].end;
    for (uint r = 1; r < records.size(); r++){
        end = end > records[r].end ? end : records[r].end;
    }

    Gene *gene        = nullptr;
    RegionType region = RegionType::NoRegion;

    if (!genome.empty()){
        if (genome.count(chrom) == 0){
            counts.outside++;
            return 1; // no genes on this chrom
        }
        if (chrom != gene_chrom){
            gene_chrom = chrom;
            set_chrom_genes(chrom, genome, genes, gene_ends);
        }

        gene = is_candidate_region(start, end, genes, gene_ends, region);
        if (gene == nullptr){
            counts.outside++;
            return 1;
        }
        if (region == RegionType::Exonic){
            counts.exonic++;
        }
        else {
            counts.promoter++;
        }
    }

    Locus locus;
    if (!build_locus(records, chrom, samples, params, locus, counts)){
        return 1;
    }

    // too much missing data in any of the 4 groups?
    if (!passes_missing(samples, locus, params.max_missing)){
        counts.too_missing++;
        return 1;
    }

    vector<uint> pure_counts[2];
    uint pure_n[2];
    tally_pure(samples, locus, pure_counts, pure_n);

    // are the two polarizing populations different (based on allele frequency)
    int ena_allele, tnp_allele;
    bool diff = classify_diff(pure_counts, pure_n, params.min_freq, ena_allele, tnp_allele);

    string pure_class = ".";
    bool shared       = false;
    if (diff){
        counts.diff++;
        pure_class = "DIFF";
        shared     = shares_homozygotes(samples, locus);
        if (shared){
            counts.shared_hom++;
            pure_class = "DIFF_SHARED_HOM";
        }
    }

    BpStats bp;
    tally_bear_paw(samples, locus, ena_allele, tnp_allele, bp);

    Classification c = Classification::Uinform;

    // are the bear paw populations heteorozygous
    if (diff && !shared && bear_paw_het_diff(bp, params.min_frac)){
        c = Classification::BpHetDiff;
        counts.bp_het_diff++;
    }
    else if (bp.het[0] + bp.het[1] > 0){
        //
        // rescue: keep any other locus with a het bear paw
        // sample so the data can still be looked at
        //
        c = Classification::BpHet;
        counts.bp_het++;
    }
    else {
        return 1; // no heterozygous bear paw samples
    }

    write_locus(ofh, locus, gene, region, c, pure_class, ena_allele, tnp_allele,
                pure_counts, pure_n, bp, samples);

    return 0;
}

int
parse_vcf(const string &vcf, const Params &params, const vector<string> &pops, 
          unordered_map<string, Sample*> &pop_map, unordered_map<string, vector<Gene*>> &genome){

    //
    // this function will be the main work horse for this code
    //
    //  1. records are read in order & grouped while their REF spans overlap
    //  2. once the next record starts past the group, the group is collapsed
    //     into a locus & tested (process_records)
    //

    bool gzipped = is_compressed(vcf);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, vcf);

    ofstream ofh("Bear_Paw_hetero_indels.tsv");

    if (!ofh.is_open()){
        cerr << "Error: Unable to create the output file in this directory\n";
        exit(EXIT_FAILURE);
    }

    ofh << std::fixed << setprecision(3);
    write_header(ofh, pops);

    //
    // sample -> pointer to Sample structs with pop & sex info
    //
    
    vector<string>  parts;
    vector<Sample*> samples;
    vector<Record>  records; // a locus can contain multiple records
    vector<Gene*>  *genes = nullptr;
    vector<uint>    gene_ends;
    string  line, chrom, rec_chrom, gene_chrom, sample;
    uint    pos, records_end = 0, last_pos = 0;
    bool    keep_reading;
    Counter counts; // singleton counter

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
        counts.records++;

        //
        // code assumes an ordered vcf
        //
        // the clustering only works on a sorted vcf
        if (chrom == rec_chrom && pos < last_pos){
            cerr << "Error: The vcf is not sorted. " << chrom << ':' << pos
                 << " comes after " << chrom << ':' << last_pos << '\n';
            exit(1);
        }
        last_pos = pos;
        Record rec;
        rec.pos = pos;
        rec.ref = parts[3];
        rec.end = pos + rec.ref.size() - 1;

        parse_string(parts[4], ',', rec.alts);

        //
        // does this record start a new locus? a new chrom or a start
        // past the furthest end of the current group means yes
        //
        if (!records.empty() && (chrom != rec_chrom || pos > records_end)){
            process_records(records, rec_chrom, genome, genes, gene_ends, gene_chrom, samples, params, ofh, counts);
            records.clear();
        }

        if (records.empty()){
            rec_chrom = chrom;
            records_end   = rec.end;
        }
        else if (rec.end > records_end){
            records_end = rec.end;
        }

        rec.parts = parts;
        records.push_back(rec);
    } // end of file parsing

    // the last group
    process_records(records, rec_chrom, genome, genes, gene_ends, gene_chrom, samples, params, ofh, counts);

    close_in_filestream(gzipped, gz_fh, txt_fh);
    ofh.close();

    cerr << "Records read " << counts.records << '\n'
         << "Loci after collapsing overlapping records " << counts.loci << '\n'
         << "  of which built from more than one record " << counts.multi_record << '\n'
         << "Loci skipped on " << MT << " & " << Y << " (haploid) " << counts.haploid_skip << '\n';

    if (!genome.empty()){
        cerr << "Loci outside the annotation regions " << counts.outside << '\n'
             << "Total exonic loci scanned "   << counts.exonic   << '\n'
             << "Total promoter loci scanned " << counts.promoter << '\n';
    }

    cerr << "Loci skipped for disagreeing REF bases " << counts.ref_mismatch << '\n'
         << "Loci skipped with only spanning deletions (*) " << counts.no_alleles << '\n'
         << "Sample calls set to missing (records could not be combined) " << counts.conflicts << '\n';

    if (params.min_gq > 0){
        cerr << "Sample calls set to missing for GQ below " << params.min_gq << ' ' << counts.low_gq << '\n';
    }
    if (params.min_ab > 0){
        cerr << "Het calls set to missing for allele balance below " << params.min_ab << ' ' << counts.low_ab << '\n';
    }

    cerr << "Loci removed for missing data " << counts.too_missing << '\n'
         << "Loci where " << pops[0] << " & " << pops[1] << " differ " << counts.diff << '\n'
         << "  of which share homozygotes " << counts.shared_hom << '\n'
         << "Loci with heterozygous bear paw for the " << pops[0] << " & " << pops[1]
         << " alleles (BP_HET_DIFF) " << counts.bp_het_diff << '\n'
         << "Loci rescued with any heterozygous bear paw sample (BP_HET) " << counts.bp_het << '\n';

    return 0;
}

void
help(){
    cerr << "Usage: ./find_bear_paw_het_indels -v indels.vcf -p popmap.tsv -m <float> -f <float> -x <float> -q <int> -b <float> -a ann.gtf [optional] -g gene_list.txt [optional] -u <int> [optional]\n"
         << "\t-v: vcf file of ONLY indels containing the bear paw samples along with the other TNP and ENA stickleback population [required]\n"
         << "\t-p: A population map containing 3 columns: [Sample][Population][Sex] [required]\n"
         << "\t-m: Minimum frequency of an allele in a population (0.5 < x <= 1), to consider it differentiated [default: 1.0]\n"
         << "\t-f: Minimum fraction of diploid Bear Paw samples in each group that must be heterozygous for the ENA & TNP alleles [default: 1.0]\n"
         << "\t-x: Maximum fraction of missing samples allowed in each group [default: 1.0; no filter]\n"
         << "\t-q: Minimum genotype quality (GQ). Calls below it are treated as missing [default: 0, off]\n"
         << "\t-b: Minimum allele balance (minor allele fraction of reads) for a het call. Hets below it are treated as missing [default: 0, off]\n"
         << "\t-a: The genome annotation file to overlap indels with\n"
         << "\t-g: A file contain the gene_ids (one per line) of interest. Only these genes will be investigated in the annotation\n"
         << "\t-u: Number of sites to look upstream of a gene when using a genome annotation [default: 0]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, pop_file, arg, genes_list;
    Params params;
    int  upstream = 0;

    // expect at least 2 arguments
    if (argc < 4){
        help();
    }

    for (int i = 1; i < argc; i++){
        arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = argv[i + 1];
            i++;
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = argv[i + 1];
            i++;
        }
        else if (arg == "-p" && i + 1 < argc){
            pop_file = argv[i + 1];
            i++;
        }
        else if (arg == "-g" && i + 1 < argc){
            genes_list = argv[i + 1];
            i++;
        }
        else if (arg == "-h"){
            help();
        }
        else if (arg == "-m" && i + 1 < argc){
            try {
                params.min_freq = stod(argv[i + 1]);

                if (params.min_freq <= 0.5 || params.min_freq > 1.0){
                    cerr << "Error: -m must be greater than 0.5 and at most 1.0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::exception &e){
                cerr << "Error: -m must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
            i++;
        }
        else if (arg == "-f" && i + 1 < argc){
            try {
                params.min_frac = stod(argv[i + 1]);

                if (params.min_frac < 0 || params.min_frac > 1.0){
                    cerr << "Error: -f must be a value within 0 and 1.0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::exception &e){
                cerr << "Error: -f must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
            i++;
        }
        else if (arg == "-x" && i + 1 < argc){
            try {
                params.max_missing = stod(argv[i + 1]);

                if (params.max_missing < 0 || params.max_missing > 1.0){
                    cerr << "Error: -x must be a value within 0 and 1.0\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::exception &e){
                cerr << "Error: -x must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
            i++;
        }
        else if (arg == "-q" && i + 1 < argc){
            try {
                params.min_gq = stoi(argv[i + 1]);

                if (params.min_gq < 0){
                    cerr << "Error: -q must be 0 or a positive integer\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::exception &e){
                cerr << "Error: -q must be an integer. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
            i++;
        }
        else if (arg == "-b" && i + 1 < argc){
            try {
                params.min_ab = stod(argv[i + 1]);

                if (params.min_ab < 0 || params.min_ab > 0.5){
                    cerr << "Error: -b must be a value within 0 and 0.5\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::exception &e){
                cerr << "Error: -b must be a float. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
            i++;
        }
        else if (arg == "-u" && i + 1 < argc){
            try {
                upstream = stoi(argv[i + 1]);

                if (upstream < 0){
                    cerr << "Error: -u must be 0 or a positive number of bp\n";
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::exception &e){
                cerr << "Error: -u must be an integer. Encountered " << argv[i + 1] << '\n';
                exit(EXIT_FAILURE);
            }
            i++;
        }
        else {
            cerr << "Unknown parameter " << arg << " provided\n";
            exit(EXIT_FAILURE);
        }
    } // end of args

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
    parse_vcf(vcf, params, pops, pop_map, genome);

    // we are done
    clean_up(genome, pop_map);

    return 0;
}