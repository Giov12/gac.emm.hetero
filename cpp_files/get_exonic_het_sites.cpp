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

using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::sort;
using std::unordered_set;
using std::unordered_map;
using std::setprecision;

//
// this code takes the output of
// either find_hwe_sites.cpp or merge_hwe.py
// & check which sites land on exons or within
// a gene
//

typedef unsigned int uint;

struct Exon {
    uint start;
    uint end;
};

struct Counts {
    uint exonic;
    uint intronic;
};

class Gene {

private:
    uint _exonic   = 0; // count how
    uint _intronic = 0; // many of each

public:
    string       id;
    uint         start;
    uint         end;
    vector<Exon> exons;

    //
    // empty constructor
    //
    Gene (string id_, uint start, uint end){
        this->id    = id_;
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

    bool has_marker(void){
        return this->_exonic > 0 || this->_intronic > 0;
    }

    void add_marker(uint marker){
        //
        // exons are ordered,
        // so just see if there is an overlap
        //

        for (uint i = 0; i < this->exons.size(); i++){
            const Exon &e = this->exons[i];
            if (e.start <= marker && marker <= e.end){
                this->_exonic++;
                return;
            }
        }

        // no overlap
        this->_intronic++;
    }

    Counts get_counts(void) const {
        //
        // return the counts of markers
        //
        Counts c;
        c.exonic   = this->_exonic;
        c.intronic = this->_intronic;

        return c;
    }
 
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
    if (gzipped){
        gz_fh = gzopen(infile.c_str(), "rb");
        if (gz_fh == NULL){
            cerr << "Error: could not open " << infile << '\n';
            exit(1);
        }
    }
    else {
        fh.open(infile);
        if (!fh.is_open()){
            cerr << "Error: could not open " << infile << '\n';
            exit(1);
        }
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

string
get_geneid(string &attributes){

    //
    // get the gene_id from a ';' delimited string
    //

    if (attributes.empty()){
        return "";
    }

    size_t start = 0, next = string::npos, length = attributes.size();
    string part;

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
                if (key == "gene_id"){
                    return value; // found it
                }
            }
       }
        // we reached the end
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }

    return "";
}

int
parse_annotation(const string &ann, unordered_map<string, vector<Gene*>> &genome, 
    const unordered_map<string, vector<uint>> &markers){
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
    Gene *g;
    string line, chrom, gene_id;
    uint start, end;
    bool eof;

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
            // grab the gene_id for this gene
            gene_id = get_geneid(parts[8]);
            if (gene_id.empty()){
                cerr << "Unable to get gene_id for the following record:\n" << line;
                exit(1);
            }
            chrom = parts[0];
            start = (uint)stoi(parts[3]);
            end   = (uint)stoi(parts[4]);
            if (parts[2] == "gene"){
                g = new Gene(gene_id, start, end);
                gene_map[chrom][gene_id] = g;
            }
            else {
                Exon exon{start, end};
                gene_map[chrom][gene_id]->add_exon(exon);
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
            chrom_genes.push_back(std::move(jtr->second));
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

int
get_markers(const string &table, unordered_map<string, vector<uint>> &markers){

    //
    // this load the markers to
    // a mapping of chrom -> <pos, pos, pos>
    //

    bool gzipped = is_compressed(table);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, table);

    vector<string> parts;
    string line, chrom;
    uint   pos, total_markers = 0, line_cnt = 0;
    bool   eof = false;

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

        line_cnt++;

        if (line_cnt == 1){
            continue; // skip header
        }

        // this should not be triggerd
        if (line.empty()){
            continue;
        }

        parse_tabular(line, parts);

        chrom = parts[0];
        pos   = (uint)stoi(parts[1]);
        markers[chrom].push_back(pos);
        total_markers++;

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);


    if (total_markers == 0){
        cerr << "Found no records in this file\n";
        return 1;
    }
    else {
        cerr << "Found a total of " << total_markers << " markers\n";
    }

    return 0;
}
void
get_overlapping_genes(const int middle, const int site, vector<Gene*> *genes, vector<Gene *> &candidates){

    //
    // populate the candidates vector with overlapping genes
    // at this site
    //
    int right = middle + 1;

    Gene *gene;

    while (right < genes->size()){
        gene = (*genes)[right];
        if (gene->start <= site && gene->end >= site){
            candidates.push_back(gene);
        }
        else if (gene->start > site){
            break;
        }
        right++;
    }
}

uint
find_overlapping(unordered_map<string, vector<uint>> &markers, unordered_map<string, vector<Gene *>> &genome){
    //
    // find the overlap between genes & the snps
    //

    vector<Gene*> *genes;
    vector<uint> *sites;
    vector<uint> gene_ends; // to handle long genes
    Gene *gene;
    string chrom;
    unordered_set<string> overlapping_genes; // should not be many

    uint overlapping = 0;

    for (auto itr = markers.begin(); itr != markers.end(); itr++){
        chrom = itr->first;

        if (genome.count(chrom) == 0){
            continue; // no genes on this chromosome
        }
        sites = &itr->second;
        genes = &genome[chrom];

        // need to handle genes with exons very spread out
        gene_ends.clear();
        gene_ends.resize(genes->size());
        uint longest = 0;
        for (uint i = 0; i < genes->size(); i++){
            longest = longest > (*genes)[i]->end ? 
                      longest : (*genes)[i]->end;
            gene_ends[i] = longest;
        }

        for (uint j = 0; j < sites->size(); j++){
            uint site = (*sites)[j];

            // binary search to find the leftmost overlapping gene
            int left  = 0, mid;
            int right = genes->size();
            while (left < right){
                mid  = left + (right - left) / 2;
                gene = (*genes)[mid];
                if (gene_ends[mid] >= site){
                    right = mid;
                }
                else {
                    left = mid + 1;
                }
            }

            bool found = false;
            if (left < genes->size()){
                // we stopped somewhere
                gene  = (*genes)[left];
                found = gene->start <= site;
                mid   = left;
            }
            if (found){
                overlapping++;

                // we need to get the range of genes this site covers
                vector<Gene *> within_range;

                within_range.push_back((*genes)[mid]);
                get_overlapping_genes(mid, site, genes, within_range);

                // now to categorize this site for each gene
                for (uint k = 0; k < within_range.size(); k++){
                    within_range[k]->add_marker(site);
                    overlapping_genes.insert(within_range[k]->id);
                }
            }
        } // end of bp loop
    } // end of chrom loop

    cerr << "Found " << overlapping << " sites overlapping " << overlapping_genes.size() << " genes\n";

    return overlapping;
}

void
write_output(unordered_map<string, vector<Gene*>> &genome){
    // write out a tsv with the genes with SNPs within their
    // boundaries

    ofstream fh("Genes_with_hetero_snps.tsv");

    fh << "#GeneID\tExonic\tIntronic\n";

    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene*> *genes = &itr->second;
        for (uint i = 0; i < genes->size(); i++){
            Gene* gene = (*genes)[i];
            if (gene->has_marker()){
                Counts c = gene->get_counts();
                fh << gene->id << '\t' << c.exonic << '\t' << c.intronic << '\n';
            }
            delete gene;
        }
    }
}

void
write_output(unordered_map<string, vector<Gene*>> &genome, const uint overlapping){
    // write out a tsv with the genes with SNPs within their
    // boundaries

    ofstream fh;

    if (overlapping > 0){
        // only write if there's something to report
        fh.open("Genes_with_hetero_snps.tsv");
        fh << "#GeneID\tExonic\tIntronic\n";
    }    

    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene*> *genes = &itr->second;
        for (uint i = 0; i < genes->size(); i++){
            Gene* gene = (*genes)[i];
            if (overlapping > 0 && gene->has_marker()){
                Counts c = gene->get_counts();
                fh << gene->id << '\t' << c.exonic << '\t' << c.intronic << '\n';
            }
            delete gene;
        }
    }

    if (overlapping > 0){
        fh.close();
    }
}

void
help(){
    cerr << "Usage: ./get_exonic_het_sites -t merged_hwe.tsv.gz -a ann.gtf.gz\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string table, ann;
    
    // expect at least two arguments
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-t" && i + 1 < argc){
            table = string(argv[i + 1]);
        }
        if (arg == "-a" && i + 1 < argc){
            ann = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (table.empty() || ann.empty()){
        help();
    }

    if (!file_exists(table)){
        cerr << "Unable to find " << table << '\n';
        exit(1);
    }
    if (!file_exists(ann)){
        cerr << "Unable to find " << ann << '\n';
        exit(1);
    }

    // first, load the markers
    unordered_map<string, vector<uint>> markers;
    get_markers(table, markers);

    // next, load genes from chromosomes with markers
    unordered_map<string, vector<Gene*>> genome;
    parse_annotation(ann, genome, markers);

    // now, find overlapping genes / markers
    uint overlapping = find_overlapping(markers, genome);

    // writes an output file only if there's something to report
    write_output(genome, overlapping);

    return 0;
}