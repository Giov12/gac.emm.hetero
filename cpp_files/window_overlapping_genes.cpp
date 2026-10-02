#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <zlib.h>

using std::string;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::unordered_map;
using std::setprecision;
using std::to_string;
using std::cerr;
using std::cout;
using std::stoi;
using std::stod;
using std::sort;

//
// code to parse the results of bootstrap_pixy
// and see which genes overlap the genomic windows
// that were significant
//

typedef unsigned int uint;

struct Gene {
    string id;
    string name;
    uint start;
    uint end;
};

struct Attribute {
    string key;
    string val;
};

struct Position {
    uint start;
    uint end;
};

struct Entry {
    string chrom;
    uint start;
    uint end;
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
parse_id(string &column, Entry &entry){

    //
    // parse a '_' delimited str
    //

    int start  = 0, end = 0, cnt = 0;
    string part;

    while (end < column.size()){
        if (column[end] == '_'){
            cnt++;
            part = column.substr(start, end - start);
            if (cnt == 1){
                entry.chrom = part;
            }
            else {
                entry.start = (uint)stoi(part);
            }
            start = end + 1;
        }
        end++;
    }

    if (start < column.size()){
        part      = column.substr(start);
        entry.end = (uint)stoi(part);
    }

    return 0;
}

void
parse_attributes(string &attributes, vector<Attribute> &atrbVec){

    //
    // parse the attributes string that is expected to
    // be ';' delimited
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
parse_annotation(const string &ann, unordered_map<string, vector<Gene>> &genome, 
                 const unordered_map<string, vector<Position>> &windows){
    //
    // collect only genes on chromosomes with markers
    //

    bool gzipped = is_compressed(ann);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, ann);

    //
    // create the objects we need to store info
    //
    vector<string> parts;
    vector<Attribute> atrbVec;
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

        if (windows.count(chrom) == 0){
            continue; // no windows on this chrom
        }

        if (parts[2] == "gene"){

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
            genome[chrom].push_back({gene_id, gene_name, start, end});
        }
    } // end of parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (genome.empty()){
        cerr << "No genes were found on chromosomes with significant window in " << ann << '\n';
        exit(1);
    }

    // move the genes into the genome map
    uint total = 0;
    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene> &genes = itr->second;
        total              += genes.size();

        // now sort for downstream binary search
        sort(genes.begin(), genes.end(),[]
            (const Gene geneA, const Gene geneB){
                if (geneA.start == geneB.start){
                    return geneA.end < geneB.end;
                }
                return geneA.start < geneB.start;
            }  
        );
    }

    cerr << "Loaded " << total << " genes\n";
    return 0;
}

int
load_windows(const string &infile, unordered_map<string, vector<Position>> &windows){
    //
    // only extract windows of significance to speed up our search
    //

    const double threshold = 0.05;

    bool gzipped = is_compressed(infile);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, infile);


    uint scanned = 0, found = 0, start, end;
    vector<string> parts;
    string line;
    double pval;
    bool eof = false;

    while (true) {
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

        if (parts.size() < 3){
            cerr << "Malformed line in " << infile << '\n' << line;
            exit(1);
        }

        if (!parts[2].empty() && parts[2].back() == '\n'){
            parts[2].pop_back();
        }
        
        scanned++;
        pval = stod(parts[2]);

        if (pval > threshold){
            continue;
        }

        // parse the id
        // chrom_start_end
        Entry entry;
        parse_id(parts[0], entry);
        windows[entry.chrom].push_back({entry.start, entry.end});
        found++;
    }

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (found == 0){
        cerr << "Did not find any windows that had a pvalue <= 0.05\n";
    }

    cerr << "Scanned " << scanned << " records. Found " << found << " windows with a p-value of <= 0.05\n";

    return 0;
}

void
merge_windows(vector<Position> &wins, vector<Position> &merged){

    //
    // sort the windows & join runs of adjacent windows that are significant
    //

    merged.clear();

    for (uint i = 0; i < wins.size(); i++){
        if (!merged.empty() && wins[i].start == merged.back().end + 1){
            merged.back().end = wins[i].end;
        }
        else {
            merged.push_back(wins[i]);
        }
    }
}

double
gene_coverage(const Gene &gene, const vector<Position> &merged, uint left, uint &covered, string &blocks){

    //
    // return the fraction of this gene that falls inside
    // significant windows & record which windows they were
    //

    uint start, end;
    blocks.clear();

    while (left < merged.size()){
        if (merged[left].start > gene.end){
            break; // windows are sorted, so we are done
        }
        start    = gene.start > merged[left].start ? gene.start : merged[left].start;
        end      = gene.end < merged[left].end ? gene.end : merged[left].end;
        covered += end - start + 1;

        if (!blocks.empty()){
            blocks += ',';
        }
        blocks += to_string(merged[left].start) + '-' + to_string(merged[left].end);
        left++;
    }

    return (double)covered / (double)(gene.end - gene.start + 1);
}


int
find_overlapping_genes(unordered_map<string, vector<Gene>> &genome, unordered_map<string, vector<Position>> &windows){
    //
    // find genes that overlap the windows
    // from the gene's perspective, 50% of the
    // gene needs to be in the window
    //

    ofstream ofh("Overlapping_genes.tsv");

    if (!ofh.is_open()){
        cerr << "Unable to create results output file in this directory\n";
        exit(1);
    }

    const double min_fraction = 0.5;
    string chrom, blocks;
    uint found = 0;

    // sort the chromosomes so the output order deterministic
    vector<string> chroms;
    for (auto itr = windows.begin(); itr != windows.end(); itr++){
        chroms.push_back(itr->first);
    }

    sort(chroms.begin(), chroms.end());

    ofh << "#Chrom\tGeneID\tGeneName\tGeneStart\tGeneEnd\tGeneLength\t"
         << "CoveredBP\tFractionCovered\tWindows\n";
    ofh << std::fixed << setprecision(3);

    for (uint i = 0; i < chroms.size(); i++){
        
        chrom = chroms[i];

        if (genome.count(chrom) == 0){
            continue; // no genes on this chrom
        }

        vector<Position> &wins = windows[chrom];
        vector<Gene>    &genes = genome[chrom];

        vector<Position> merged; // merge overlapping significant windows
        merge_windows(wins, merged);

        for (uint g = 0; g < genes.size(); g++){
            const Gene &gene = genes[g];

            // binary search for the first merged window that
            // ends at or after the gene's start
            uint left = 0, right = merged.size();
            while (left < right){
                uint mid = left + (right - left) / 2;
                if (merged[mid].end >= gene.start){
                    right = mid;
                }
                else {
                    left = mid + 1;
                }
            }

            uint covered = 0;
            double fraction = gene_coverage(gene, merged, left, covered, blocks);

            if (fraction < min_fraction){
                continue;
            }
            uint length = gene.end - gene.start + 1;
            fraction    = fraction * 100.0;

            ofh  << chrom     << '\t' << gene.id    << '\t'
                 << gene.name << '\t' << gene.start << '\t' 
                 << gene.end  << '\t' << length     << '\t'
                 << covered   << '\t' << fraction   << '\t' 
                 << blocks    << '\n';
            found++;
        }
    }

    cerr << "Reported " << found << " genes with >= "
         << (min_fraction * 100) << "% of their length in significant windows\n";


    ofh.close();

    return 0;
}

void
help(void){
    cerr << "Usage: ./window_overlapping_genes -f Windows_pvalues.tsv -a an.gtf\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string infile, ann;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-f" && i + 1 < argc){
            infile = argv[i + 1];
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = argv[i + 1];
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (infile.empty() || ann.empty()){
        help();
    }

    if (!file_exists(infile)){
        cerr << "Unable to find " << infile << '\n';
        exit(1);
    }

    if (!file_exists(ann)){
        cerr << "Unable to find " << ann << '\n';
        exit(1);
    }

    //
    // get signficant windows first
    //
    unordered_map<string, vector<Position>> windows;
    load_windows(infile, windows);
   
    //
    // now get the genes
    //
    unordered_map<string, vector<Gene>> genome;
    parse_annotation(ann, genome, windows);

    // now overlap the structures
    find_overlapping_genes(genome, windows);
    
    return 0;
}