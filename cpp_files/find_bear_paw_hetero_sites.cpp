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
using std::unordered_map;
using std::unordered_set;
using std::setprecision;
using std::stod;

//
// this code is specifically designed. for finding heterozygous
// sites for the population of interest (bear paw). If given an
// annotation file, only exonic sides will be searched. If given a
// list of genes, then only those genes within the annotation will be
// searched
//

typedef unsigned int uint;

enum Sex {
    Male   = 0,
    Female = 1
};

enum mtType {
    ENA = 0,
    TNP = 1
};

enum Classification {
    Uinform = 0, // un-informative
    Diff    = 1,
    Private = 2,
};
struct Exon {
    uint start;
    uint end;
};
struct Attribute {
    string key;
    string val;
};

class Gene {
private:
    bool _coding = false;

public:
    string       id;
    string       name;
    uint         start;
    uint         end;
    vector<Exon> exons;

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

    bool is_exonic(uint pos){
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

    void set_as_coding(void){
        this->_coding = true;
    }

    bool is_coding(void) const {
        return this->_coding;
    }
 
};

struct Sample {
    mtType pop;
    Sex  sex;
    bool bear_paw;
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
                 const bool only_targets, unordered_set<string> &target_genes){
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
                g = new Gene(gene_id, gene_name, start, end);
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
            chrom_genes.push_back(jtr->second);
            total++;
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
is_exonic(const uint pos, vector<Gene*>*genes, vector<uint> &gene_ends){
    //
    // helper function to determine
    // whether this site lands on a single exon
    // return a nullptr if no overlap
    //

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
    bool found = false; // did we find a candidate
    if (left < ngenes){
        gene  = (*genes)[left];
        found = gene->start <= pos;
    }

    if (found){
        while (left < ngenes){
            gene = (*genes)[left];
            if (gene->start <= pos && pos <= gene->end){
                if (gene->is_exonic(pos)){
                    return gene; // lands in an exon
                }
            }
            else if (pos < gene->start){
                break; // position is before the gene begins
            }
            left++;
        }
    }

    return nullptr; // never hit true
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
is_informative(const uint bp[2][4], const double min_frac){

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

    double frac_ena = (double)(bp[0][0] + bp[0][1]) / called_ena;
    double frac_tnp = (double)(bp[1][2] + bp[1][1]) / called_tnp;

    return (frac_ena >= min_frac) && (frac_tnp >= min_frac);
}

void
write_header(ofstream &fh, const vector<string> &pops){
    //
    // write the column headers for the output tsv
    //
    fh << "#Chrom\tBP\tGene_id\tRef\tAlt\tClass\t"
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
        summary += std::to_string(counts[i]) + ' ' + labels[i];
    }

    if (summary.empty()){
        summary = "."; // no samples?
    }

    return summary;
}

void
write_site(ofstream &fh, const string &chrom, const uint pos, const Gene *gene,
           const string &ref, const string &alt, const Classification c,
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

    string c_label;
    if (c == Classification::Diff){
        c_label = "DIFF";
    }
    else if (c == Classification::Private){
        c_label = "PRIVATE";
    }
    else {
        c_label = "UNINFORMATIVE"; // should not happen
    }

    string gene_id = ".";
    if (gene != nullptr){
        gene_id = gene->id;
    }

    fh << chrom    << '\t' << pos      << '\t' << gene_id << '\t'
       << ref      << '\t' << alt      << '\t' << c_label << '\t'
       << ena_base << '\t' << tnp_base << '\t'
       << ena_freq << '\t' << tnp_freq << '\t'
       << pure_n[0] << '\t' << pure_n[1];

    for (uint i = 0; i < 2; i++){
        for (uint j = 0; j < 4; j++){
            fh << '\t' << bp[i][j];
        }
    }
    // genotype summaries: ENA, BP ENA, BP TNP, TNP
    fh << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::ENA, false)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::ENA, true)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::TNP, true)
       << '\t' << genotype_summary(samples, ploidy, alt_copies, ref, alt, mtType::TNP, false) << '\n';

}

int
parse_vcf(const string &vcf, const bool only_exonic, const double min_freq,
          const double min_frac, const double max_missing, const vector<string> &pops,
          unordered_map<string, Sample*> &pop_map, unordered_map<string, vector<Gene*>> &genome){

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
    
    //
    // sample -> pointer to Sample structs with pop & sex info
    //
    
    vector<string>  parts;
    vector<Gene*>  *genes;
    vector<Sample*> samples;
    vector<uint>    gene_ends, ploidy;
    vector<int>     alt_copies;
    Classification clsfion;
    string line, sample, chrom, curChrom;
    uint   pos, ngenes, furthest;
    bool   keep_reading;
    int    gt_idx, ena_allele;
    uint valid_sites = 0, snp_sites = 0, informative = 0, passed = 0, miss_cnt = 0;
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

        chrom = parts[0];
        pos   = (uint)stoi(parts[1]);
        gene  = nullptr;

        // if using only genes, confirm that it is exonic about the exonic regions
        if (only_exonic){ 

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
                    if ((*genes)[n]->end > furthest){
                        furthest = (*genes)[n]->end; 
                    }
                    gene_ends[n] = furthest;
                }
            }
            
            // check if this position is exonic
            gene = is_exonic(pos, genes, gene_ends);
            if (gene == nullptr){
                continue;
            }
            valid_sites++;
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

        if (ena_allele == -1){
            continue;;
        }

        informative++;

        // does the bear paw samples carry their mtDNA's type?
        tally_bear_paw(samples, ploidy, alt_copies, ena_allele, bp);

        if (!is_informative(bp, min_frac)){
            continue;
        }
        passed++;
        write_site(ofh, chrom, pos, gene, ref, alt, clsfion, ena_allele, pure_alt, pure_n, bp, 
                   samples, ploidy, alt_copies);


    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);
    ofh.close();

    if (only_exonic){
        cerr << "Total exonic sites scanned " << valid_sites << '\n';
    }

    cerr << "Biallelic SNPs scanned "    << snp_sites << '\n'
         << "Sites removed for missing data " << miss_cnt << '\n'
         << "Informative sites between " << pops[0]   << " and " << pops[1] << ' ' << informative << '\n'
         << "Sites passing the bear paw check " << passed << '\n';

    return 0;
}

void
help(){
    cerr << "Usage: ./find_bear_paw_hetero_sites -v vcf_file -p popmap.tsv -m <float> -f <float> -x <float> -a ann.gtf [optional] -g gene_list.txt [optional]\n"
         << "\t-v: vcf file containing the bear paw samples along with the other TNP and ENA stickleback population [required]\n"
         << "\t-p: A population map containing 3 columns: [Sample][Population][Sex] [required]\n"
         << "\t-m: Minimum frequency of an allele in a population (0.5 < x <= 1), to consider it differentiated [default: 1.0]\n"
         << "\t-f: Minimum fraction of Bear Paw Samples to carry the specific mt-Type allele [default: 1.0]\n"
         << "\t-x: Maximum fraction of missing samples allowed in each group [default: 1.0; no filter]\n"
         << "\t-a: The genome annotation file to overlap SNPs with\n"
         << "\t-g: A file contain the gene_ids (one per line) of interest. Only these genes will be investigated in the annotation\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, pop_file, arg, genes_list;
    double min_freq = 1.0, min_frac = 1.0, max_missing = 1.0;
    
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
        parse_annotation(ann, genome, only_targets, target_genes);
    }

    // load the samples
    unordered_map<string, Sample*> pop_map;
    vector<string> pops; // index 0 -> ENA, index 1 -> TNP
    parse_popmap(pop_file, pop_map, pops);
    
    // parse the vcf file
    parse_vcf(vcf, use_ann, min_freq, min_frac, max_missing, pops, pop_map, genome);

    // we are done
    clean_up(genome, pop_map);

    return 0;
}