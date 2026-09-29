#include <iostream>
#include <cstring>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <math.h>
#include <sys/stat.h>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::stoi;
using std::stod;
using std::ifstream;
using std::ofstream;

//
// this code will just make any outlier Fst 
// regions standout relative to the
// genomic background
//

typedef unsigned int uint;

// sigma & weights will be global
const int SIGMA            = 150000;
const int window_span      = (3 * SIGMA) + 1;
const int buffer_line_size = 1024; 
double weights[window_span];

struct Chrom {
    string         name;
    vector<uint>   positions;
    vector<double> fsts;
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
calc_weights(double (&weights)[window_span]){

    // fill weights once
    for (int i = 0; i < window_span; i++){
        weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(SIGMA, 2)));
    }
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
smooth_fsts(Chrom &current_chrom, ofstream &filestream){
    
    /*
    positions      = bp positions
    fsts           = fsts at each position
    weights        = precomputed kernel smoothed weights
    start_pos      = beginning of window on chrom
    end_pos        = end of window on chrom
    window_cen     = current position on chrom
    window_start   = used to find start_pos by adjusting from center
    window_end     = used to find end_pos by adjusting to window_start
    n_window_sites = number of positions
    distnace       = distance from window center
    window_avg     = will hold weighted average
    window_weight  = will accumulate weights to use for weighing average
    weight         = kernel weight at position
    */

    int start_pos = 0, end_pos = 0;
    const int n_sites = current_chrom.positions.size();
    int window_cen, window_start, window_end, n_window_sites, distance, pos_num;
    double avg, window_weight, weight, fst;
    vector<double> smoothed_fsts(n_sites, 0.0);

    if (n_sites == 0){
        return;
    }

    for (int site = 0; site < n_sites; site++){
        window_cen   = current_chrom.positions[site];       // go to center
        window_start = window_cen - (window_span / 2);      // go to start
        window_start = window_start < 1 ? 1 : window_start; // keep within bounds
        window_end   = window_cen + (window_span / 2) + 1;
        while (current_chrom.positions[start_pos] < window_start){ // start pointer in positions
            start_pos++;
        }
        while (end_pos < n_sites && current_chrom.positions[end_pos] < window_end){ // end pointer in positions
            end_pos++;
        }
        n_window_sites = 0; // window size
        avg            = 0.0; // reset average
        window_weight = 0.0;
        for (int pos = start_pos; pos < end_pos; pos++){
            n_window_sites++;
            fst      = current_chrom.fsts[pos];
            pos_num  = current_chrom.positions[pos]; // get current position
            distance = (int)abs(pos_num - window_cen); // relative to center
            weight   = weights[distance]; // get weight at position
            avg     += (fst * weight); // accumulate sum
            window_weight += weight; // accumulate weight
        }
        if (n_window_sites == 0){
            exit(1);
        }
        smoothed_fsts[site] = avg/window_weight;
    }

    // now to write to output
    for (int i = 0; i < n_sites; i++){
        filestream << current_chrom.name << '\t' << current_chrom.positions[i] << '\t' 
                   << smoothed_fsts[i] << '\n';
    }

    // clear chrom struct for next chromosome
    current_chrom.name.clear();
    current_chrom.positions.clear();
    current_chrom.fsts.clear();
}

string 
make_output_name(const string &fst_table){

    // simple helper function to create an output name

    // get basename of file
    string file_basename = fst_table.substr(fst_table.find_last_of("/\\") + 1);

    // add suffix
    return "Smoothed_" + file_basename.substr(0, file_basename.find_last_of('.'));
}

int
smooth_fst_file(string &fst_file){

    //
    // find snps that are overlapping genes & add each populations
    // score to each gene
    //

    bool gzipped = is_compressed(fst_file);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    ofstream ofh(make_output_name(fst_file));
    open_in_filestream(gzipped, gz_fh, txt_fh, fst_file);

    vector<string> parts;
    string line, chrom;
    Chrom curChrom; // current chromosome
    double fst;
    uint   pos;
    bool eof = false;

    //
    // skip the header
    //
    uint line_num = 0;

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
        line_num++;

        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }

        parse_tabular(line, parts);

        if (line_num == 1 && parts.front() == "CHROM"){
            continue; // skip header
        }

        if (parts.size() < 3){
            cerr << "Error: Malformed line found:\n" << line << '\n';
            exit(1);
        }

        chrom = parts[0];
        if (curChrom.name.empty()){
            curChrom.name = chrom;
        }
        if (parts[2][0] == '-'){
            continue; // these are effectively 0
        }
        fst = stod(parts[2]);
        pos = (uint)stoi(parts[1]);

        if (curChrom.name == chrom){
            curChrom.fsts.push_back(fst);
            curChrom.positions.push_back(pos);
        }
        else {
            smooth_fsts(curChrom, ofh);
            curChrom.name = chrom; // reset
            curChrom.fsts.push_back(fst);
            curChrom.positions.push_back(pos);
        }

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    // now process the last chromosomes
    if (!curChrom.fsts.empty()){
        smooth_fsts(curChrom, ofh);
    }

    ofh.close();

    return 0;
}

void
help(void){
    cerr << "./kernel_smooth_vcftools_fst -f fst_file [required]\n";
    exit(1);
}
int main(int argc, char *args[]){

    // params
    string fst_file, param;

    for (int i = 0; i < argc; i++){
        param = string(args[i]);

        if (param == "-h" || param == "--help") {
            help();
        }
        else if (param == "-f" && i + 1 < argc){
            fst_file = string(args[i + 1]);
        }
    }

    if (fst_file.empty()){
        cout << "No fst file supplied.\n";
        exit(1);
    }

    // compute the weights
    calc_weights(weights);

    // s single function will do the heavy lifting
    smooth_fst_file(fst_file);
    
    return 0;
}