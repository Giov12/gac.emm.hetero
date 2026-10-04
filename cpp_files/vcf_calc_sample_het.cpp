#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <zlib.h>
#include <iomanip>
#include <algorithm>
#include <unordered_map>

using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::sort;
using std::unordered_map;
using std::setprecision;

//
// some code I wrote just to get an idea on the amount
// heterozygous sites per sample with an optional
// feature to only calculate heterozygosity in exonic
// sites at variant positions
//

typedef unsigned int uint;

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
    string sample;
    uint   tot; // total sites
    uint   het; // hetero sites
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
parse_annotation(const string &ann, unordered_map<string, vector<Gene*>> &genome){
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
                continue; // skip non-coding genes
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
is_hetero(string &geno){
    
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

    return a1 == a2 ? 0 : 1;
}

bool
is_exonic(const uint pos, vector<Gene*>*genes, vector<uint> &gene_ends){
    //
    // helper function to determine
    // whether this site lands on a single exon
    //

    if (genes->empty()){
        return false; // no genes to search through
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
                    return true; // lands in an exon
                }
            }
            else if (pos < gene->start){
                break; // position is before the gene begins
            }
            left++;
        }
    }

    return false; // never hit true
}

int
parse_vcf(const string &vcf, const bool only_genic,
          unordered_map<string, vector<Gene*>> &genome){

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
    vector<Sample> samples;
    vector<Gene*> *genes;
    vector<uint> gene_ends;
    string line, sample, chrom, curChrom;
    uint   pos, ngenes, furthest;
    bool   eof = false; // default value

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
                    samples.push_back({sample, 0, 0});
                }
            }
            continue;
        }

        if (samples.empty()){
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
        parse_tabular(line, parts);

        if (only_genic){ // we care only about the exonic regions
            chrom = parts[0];

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
            pos = (uint)stoi(parts[1]);

            // this is not an exonic site
            if (!is_exonic(pos, genes, gene_ends)){
                continue;
            }
        }

        for (uint i = 9; i < parts.size(); i++){
            uint sample_idx = i - 9;
            string geno     = parts[i];
            geno            = geno.substr(0, geno.find(':')); 
            if (geno.find('.') != string::npos){ // missing data
                continue;
            }
            samples[sample_idx].tot++;
            samples[sample_idx].het += is_hetero(geno);
        }

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    cerr << std::fixed << std::setprecision(2);
    
    for (uint i = 0; i < samples.size(); i++){
        Sample &s = samples[i];
        double hetero = s.tot > 0 ? ((double)s.het / (double)s.tot) * 100.0 : 0.0;
        cerr << "Sample " << s.sample << ' ' << hetero << "% hetero\n";
    }

    return 0;
}

void
help(){
    cerr << "Usage: ./vcf_calc_sample_het -v vcf_file -a ann.gtf [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, arg;
    unordered_map<string, vector<Gene*>> genome;
    
    // expect at least a single argument
    if (argc < 2){
        help();
    }

    for (int i = 1; i < argc; i++){
        arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = argv[i + 1];
        }
        if (arg == "-a" && i + 1 < argc){
            ann = argv[i + 1];
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (vcf.empty()){
        help();
    }

    if (!file_exists(vcf)){
        cerr << "Unable to find " << vcf << '\n';
        exit(1);
    }

    //
    // check if we will use an annotation file
    // to only look at exonic regions for
    // protein coding genes
    //
    bool use_ann = !ann.empty();
    if (use_ann){
        if (!file_exists(ann)){
            cerr << "Unable to find " << ann << '\n';
            exit(1);
        }
        parse_annotation(ann, genome);
    }
    
    // parse the vcf file
    parse_vcf(vcf, use_ann, genome);

    return 0;
}