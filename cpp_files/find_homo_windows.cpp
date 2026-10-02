#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <charconv>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <zlib.h>

using std::string;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;

//
// code to parse the results of merge_hwe.py
// and find regions of the genome where the
// population of interest is homozygous
//

typedef unsigned int uint;

struct Entry {
    int homo_ref = -1; // may not be necessary but as a safe guard
    int hetero   = -1;
    int homo_alt = -1;
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
parse_column(string &column, Entry &e){

    //
    // parse a '/' delimited line
    //

    int start  = 0, end = 0, cnt = 0, val;

    while (end < column.size()){
        if (column[end] == '/'){
            val = stoi(column.substr(start, end - start));
            cnt++;
            
            switch (cnt)
                {
                case 1:
                    e.homo_ref = val;
                    break;
                case 2:
                    e.hetero = val;
                    break;
                case 3:
                    e.homo_alt = val; // should not happen
                    break;
                default:
                    break;
                }
            start = end + 1;
        }
        end++;
    }

    if (start < column.size()){
        val = stoi(column.substr(start, end - start));
        if (cnt != 2){
            cerr << "Bad column found: " << column << '\n';
            exit(1);
        }
        e.homo_alt = val;
    }

    return 0;
}

bool
is_homozygous(vector<string> &columns, const int index){
    //
    // is this site a candidate region
    //

    // first check if the focal population is homozygous
    Entry pop_entry;
    parse_column(columns[index], pop_entry);

    if (pop_entry.hetero > 0){
        return false; // heterozygous
    }
    else if (pop_entry.homo_alt == 0 && pop_entry.homo_ref == 0){
        return false; // site is absent of data
    }
    else if (pop_entry.homo_alt > 0 && pop_entry.homo_ref > 0){
        return false; // homozygous for both alleles
    }

    // figure out which allele this pop is homozgyous for
    bool ref_homo = pop_entry.homo_ref > 0;
    bool alt_homo = pop_entry.homo_alt > 0;

    // now check the other columns
    for (int i = 2; i < columns.size(); i++){
        if (i == index){
            continue; // skip self
        }
        Entry e;
        parse_column(columns[i], e);
        if (ref_homo){
            if (e.homo_ref > 0){
                return false;
            }
        }
        else if (alt_homo){
            if (e.homo_alt > 0){
                return false;
            }
        }
    }

    return true; // could not disprove it

}

int
slide_window(const string &chrom, vector<int> &positions, const int window, ofstream &ofh){

    //
    // helper function to slide through the chromosome and find windows of interest
    // function assumes positions is already ordered
    //

    if (positions.empty()){
        return 1; // edge case
    }

    int start, prev;

    for (uint i = 0; i < positions.size(); i++){
        if (i == 0){
            start = positions[i]; // set as default
            prev  = positions[i];
            continue;
        }
        else if (prev + window >= positions[i]){
            prev = positions[i];
        }
        else { 
            if (prev - start + 1 > 1){
                ofh << chrom << '\t' << start << '\t' << prev << '\t' << (prev - start + 1) << '\n';
            }
            start = positions[i];
            prev  = positions[i];
        }
    }

    if (prev - start + 1 > 1){
        ofh << chrom << '\t' << start << '\t' << prev << '\t' << (prev - start + 1) << '\n';
    }

    positions.clear();

    return 0;
}

int
parse_table(const string &infile, const string &pop, const int window){
    //
    // only extract windows of significance to speed up our search
    //

    const double threshold = 0.05;

    bool gzipped = is_compressed(infile);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, infile);

    ofstream ofh("Heterozygous_windows.tsv");

    if (!ofh.is_open()){
        cerr << "Unable to create output file in this directory\n";
        exit(1);
    }

    vector<int> positions;
    vector<string> parts;
    int pop_index = -1, pos;
    string line, curChrom, chrom;
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

        if (line.empty()){
            continue;
        }

        if (line.back() == '\n'){
            line.pop_back();
        }

        parse_tabular(line, parts);

        if (parts.size() < 3){
            cerr << "Malformed line in " << infile << '\n' << line;
            exit(1);
        }

        if (pop_index == -1){
            // we need to get the population of interest
            for (uint i = 2; i < parts.size(); i++){
                if (parts[i] == pop){
                    pop_index = i;
                    break;
                }
            }
            if (pop_index == -1){
                cerr << "Unable to find " << pop << " in header\n";
                exit(1);
            }
            continue;
        }


        chrom = parts[0];
        pos   = stoi(parts[1]);

        if (curChrom != chrom){
            // flush out & update
            slide_window(chrom, positions, window, ofh);
            curChrom = chrom;
        }

        if (is_homozygous(parts, pop_index)){
            positions.push_back(pos);
        }
    }

    // get last bit
    if (!positions.empty()){
        slide_window(chrom, positions, window, ofh);
    }

    close_in_filestream(gzipped, gz_fh, txt_fh);
    ofh.close();

    return 0;
}

int
create_int(const char *arg){

    // if possible, create a number out of this
    // char array


    if (arg == nullptr || *arg == '\0'){
        cerr << "Empty input provided for -w [sliding window size] ";
        exit(1);
    }

    const char *end = arg + std::strlen(arg);
    int val         = 0;
    auto result     = std::from_chars(arg, end, val, 10);

    if (result.ec != std::errc() || result.ptr != end){
        cerr << arg << " is not a valid sliding window [-w] value\n";
        exit(1);
    }

    return val;
}

void
help(void){
    cerr << "Usage: ./find_homo_windows -t merged_hwe.tsv.gz -p POP -w INT [default 250]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string infile, arg, pop;
    int window = 250;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        arg = argv[i];
        if (arg == "-t" && i + 1 < argc){
            infile = argv[i + 1];
        }
        else if (arg == "-p" && i + 1 < argc){
            pop = argv[i + 1];
        }
        else if (arg == "-w" && i + 1 < argc){
            window = create_int(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (infile.empty()){
        help();
    }

    if (!file_exists(infile)){
        cerr << "Unable to find " << infile << '\n';
        exit(1);
    }

    //
    // a single function will orchestrate everything
    //
    parse_table(infile, pop, window);
    
    return 0;
}