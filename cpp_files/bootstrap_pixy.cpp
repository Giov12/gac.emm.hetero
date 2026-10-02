#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <charconv>
#include <cstring>
#include <cstdint>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <zlib.h>
#include <random>
#include <omp.h> // requires -fopenmp

using std::string;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::stod;
using std::sort;
using std::setprecision;
using std::seed_seq;
using std::mt19937_64;
using std::uniform_int_distribution;
using std::lower_bound;

//
// code to parse the results of a table generated from pixy that will
// bootstrap the scores and return a p-value
// of whether the value is significant based on how
// often it was sampled across replicates
//

typedef unsigned int uint;

struct Window {
    string id;
    double score;
    double pval;
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
parse_pixy(string &table, string &pop1, string &pop2, vector<Window> &windows, const bool skip){

    //
    // parse the results of a pixy output file and load the windows
    // for the population comparison of interest
    //

    bool gzipped = is_compressed(table);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, table);


    vector<string> parts;
    string line, chrom;
    bool eof = false;

    //
    // counters to find the column &
    // number of variant sites within exons
    //
    uint found = 0, line_num = 0;

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

        if (line_num == 1){
            if (parts[0] != "pop1"){
                cerr << "Did not encounter a header starting with pop1 in " << table << '\n';
                exit(1);
            }
            continue; // skip header
        }

        if (parts.size() < 7){
            cerr << "Error: Expected at least 7 columns. Offending line: " << line << '\n';
            exit(1);
        }

        // check that this is for the populations of interest
        bool valid = false; 

        if (parts[0] == pop1 && parts[1] == pop2){
            valid = true;
        }
        else if (parts[0] == pop2 && parts[1] == pop1){
            valid = true;
        }

        if (!valid){
            continue;
        }

        chrom = parts[2];
        Window window;
        window.id = chrom + '_' + parts[3] + '_' + parts[4];

        if (parts[5] == "NA" || parts[5][0] == '-'){
            if (skip && parts[5] == "NA"){
                continue;
            }
            window.score = 0.0;
        }
        else {
            window.score = stod(parts[5]);
        }
      
        windows.push_back(window);
    }

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (windows.empty()){
        cerr << "Did not find any windows between " << pop1 << " and " << pop2 << '\n';
        exit(1);
    }

    cerr << "Loaded " << windows.size() << " windows between " << pop1 << " and " << pop2 << '\n';

    return 0; 
}
int
bootstrap(vector<Window> &windows, const uint bootstraps, const uint threads, const uint seed){

    //
    // function that will actually compute the bootstrap replicates in parallel
    //

    const uint nwindows = windows.size();
    uint rounds       = 0;

    // sort the values
    vector<double> sorted_scores(nwindows);

    for (uint i = 0; i < nwindows; i++){
        sorted_scores[i] = windows[i].score;
    }

    sort(sorted_scores.begin(), sorted_scores.end());

    //
    // how often a score in sorted_scores[i] was seen
    // across replicates
    //
    vector<uint> counts(nwindows, 0);

    #pragma omp parallel num_threads(threads)
    {

        vector<uint> local_counts(nwindows, 0);

        #pragma omp for schedule(static)
        for (uint b = 0; b < bootstraps; b++){
            // set the seed
            seed_seq ss{seed, b};

            // create the random generator
            mt19937_64 rng(ss);

            // create a random sampler by sampling the positions of each score
            uniform_int_distribution<uint> pick(0, nwindows - 1);

            // increment the counts
            for (uint p = 0; p < nwindows; p++){
                local_counts[pick(rng)]++;
            }

            uint finished;
            #pragma omp atomic capture
            finished = ++rounds;

            if (finished % 1000 == 0){
                #pragma omp critical(log)
                cerr << "Finished round " << finished << '\n';
            }
        }
        #pragma omp critical
        {
            for (uint j = 0; j < nwindows; j++){
                counts[j] += local_counts[j];
            }
        }
    } // end of pragma

    //
    // start counting backwards
    // idea: how many score values drawn at random
    // are greater than the score value at position i
    //
    vector<uint64_t> ranks(nwindows + 1, 0);
    for (uint i = nwindows; i-- > 0;){ // ensure we do not hit nwindows + 1
        ranks[i] = ranks[i + 1] + counts[i];
    }

    const double total = (double)bootstraps * (double)nwindows;

    // now assign the p-values
    for (uint i = 0; i < nwindows; i++){ // binary search to find the first position (iterator) that is >= than windows[i].score
        uint j          = lower_bound(sorted_scores.begin(), sorted_scores.end(), windows[i].score) - sorted_scores.begin();
        windows[i].pval = ((ranks[j] + 1.0)/( total + 1.0)); // add 1.0 to prevent things ever being zero
    }

    return 0;
}

int
write_output(vector<Window> &windows){
 
    //
    // write a simple 3-column tsv
    // where window ID, score, pvalue
    //

    string outname = "Window_pvalues.tsv";
    ofstream fh(outname);

    if (!fh.is_open()){
        cerr << "Error: Unable to write " << outname << " in this directory\n";
        exit(1);
    }

    fh << "#Window\tValue\tpvalue\n";
    fh << setprecision(6);

    for (uint i = 0; i < windows.size(); i++){
        fh << windows[i].id << '\t' << windows[i].score << '\t' << windows[i].pval << '\n';
    }

    fh.close();

    return 0;
}

uint
create_uint(const char *arg, const uint n){

    // if possible, create a number out of this
    // char array

    string param;
    switch (n)
    {
    case 0:
        param = "--bootstraps";
        break;
    case 1:
        param = "--threads";
        break;
    case 2:
        param = "--seed";
        break;
    default:
        break;
    }

    if (arg == nullptr || *arg == '\0'){
        cerr << "Empty input provided for " << param << '\n';
        exit(1);
    }

    const char *end = arg + std::strlen(arg);
    uint val        = 0;
    auto result     = std::from_chars(arg, end, val, 10);

    if (result.ec != std::errc() || result.ptr != end){
        cerr << arg << " is not a valid " << param << " value\n";
        exit(1);
    }

    return val;
}

void
help(){
    cerr << "Usage: ./bootstrap_pixy -f pixy_results.txt --pop1 POP1 --pop2 POP2 --bootstraps INT [default: 10000] --threads INT [default 1] --seed INT [default 1234] --skip_nodata [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string infile, pop1, pop2;
    uint bootstraps  = 10000;
    uint threads     = 1;
    uint seed        = 1234;
    bool skip_nodata = false;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-f" && i + 1 < argc){
            infile = argv[i + 1];
        }
        else if (arg == "--skip_nodata"){
            skip_nodata = true;
        }
        else if (arg == "--bootstraps" && i + 1 < argc){
            bootstraps = create_uint(argv[i + 1], 0);
        }
        else if (arg == "--threads" && i + 1 < argc){
            threads = create_uint(argv[i + 1], 1);
        }
        else if (arg == "--seed" && i + 1 < argc){
            seed = create_uint(argv[i + 1], 2);
        }
        else if (arg == "--pop1" && i + 1 < argc){
            pop1 = string(argv[i + 1]);
        }
        else if (arg == "--pop2" && i + 1 < argc){
            pop2 = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (infile.empty()){
        help();
    }

    if (pop1.empty()){
        cerr << "--pop1 is required\n";
        exit(1);
    }

    if (pop2.empty()){
        cerr << "--pop2 is required\n";
        exit(1);
    }

    if (!file_exists(infile)){
        cerr << "Unable to find " << infile << '\n';
        exit(1);
    }

    if (bootstraps < 1){
        cerr << "Bootstraps must be at least 1\n";
        exit(1);
    }

    if (threads < 1){
        cerr << "Threads must be at least 1\n";
        exit(1);
    }
    else if (threads > (uint)omp_get_max_threads()){
        cerr << "Warning: Max threads available is " << omp_get_max_threads() << '\n';
        threads = (uint)omp_get_max_threads();

    }

    // read in the data
    vector<Window> windows;
    parse_pixy(infile, pop1, pop2, windows, skip_nodata);

    cerr << "Starting " << bootstraps << " bootstraps\n";

    // next, count how often we see each score
    // value and rank it's frequency based on
    // how many times we see other scores higher than it
    bootstrap(windows, bootstraps, threads, seed);

    // write the output
    write_output(windows);
    
    return 0;
}