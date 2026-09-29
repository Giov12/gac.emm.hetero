#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <charconv>
#include <cstring>
#include <algorithm>
#include <vector>
#include <zlib.h>
#include <random>
#include <omp.h>

using std::string;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::stof;
using std::stod;
using std::sort;
using std::setprecision;
using std::seed_seq;
using std::mt19937_64;
using std::uniform_int_distribution;
using std::lower_bound;


//
// code to parse the results of calc_gene_fst that will
// bootstrap the Fst scores and return a p-value
// of whether the value is significant based on how
// often it was sampled across replicates
//

typedef unsigned int uint;

struct Gene {
    string id;
    float fst;
    float pval;
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
load_genes(const string &infile, vector<Gene> &genes, const bool transcript_level){

    //
    // find snps that are overlapping genes & add each populations
    // score to each gene
    //

    bool gzipped = is_compressed(infile);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, infile);


    vector<string> parts;
    string line, gene;
    float fst;
    const uint column = transcript_level ? 2 : 1;
    bool eof      = false;
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

        line_num++;

        if (line_num == 1 && !line.empty() && line[0] == '#'){
            continue; // skip header
        }


        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }

        parse_tabular(line, parts);

        if (parts.size() != 5){
            cerr << "Error: Malformed line found:\n" << line << '\n';
            exit(1);
        }
        gene = parts[column];
        fst  = stof(parts[4]);
        genes.push_back({gene, fst, 1.0}); // set p-value to 1.0 as default

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (genes.empty()){
        cerr << "No genes loaded from " << infile << '\n';
        exit(1);
    }

    cerr << "Loaded " << genes.size() << " genes\n";

    return 0;
}

int
boostrap(vector<Gene> &genes, const uint boostraps, const uint threads, const uint seed){

    //
    // function that will actually compute the bootstrap replicates in parallel
    //

    const uint ngenes = genes.size();

    // sort the fst values
    vector<float> sorted_fsts(ngenes);

    for (uint i = 0; i < ngenes; i++){
        sorted_fsts[i] = genes[i].fst;
    }

    sort(sorted_fsts.begin(), sorted_fsts.end());

    //
    // how often a fst in sorted_fsts[i] was seen
    // across replicates
    //
    vector<uint> counts(ngenes, 0);

    #pragma omp parallel num_threads(threads)
    {

        vector<uint> local_counts(ngenes, 0);

        #pragma omp for schedule(static)
        for (uint b = 0; b < boostraps; b++){
            // set the seed
            seed_seq ss{seed, b};

            // create the random generator
            mt19937_64 rng(ss);

            // create a random sampler by sampling the positions of each fst
            uniform_int_distribution<uint> pick(0, ngenes - 1);

            // increment the counts
            for (uint p = 0; p < ngenes; p++){
                local_counts[pick(rng)]++;
            }
        }
        #pragma omp critical
        {
            for (uint j = 0; j < ngenes; j++){
                counts[j] += local_counts[j];
            }
        }
    } // end of pragma

    //
    // start counting backwards
    // idea: how many fst values drawn at random
    // are greater than the fst value at position i
    //
    vector<uint64_t> ranks(ngenes + 1, 0);
    for (uint i = ngenes; i-- > 0;){ // ensure we do not hit ngenes + 1
        ranks[i] = ranks[i + 1] + counts[i];
    }

    const double total = (double)boostraps * (double)ngenes;

    // now assign the p-values
    for (uint i = 0; i < ngenes; i++){ // binary search to find the first position (iterator) that is > than genes[i].fst
        uint j        = lower_bound(sorted_fsts.begin(), sorted_fsts.end(), genes[i].fst) - sorted_fsts.begin();
        genes[i].pval = (float)((ranks[j] + 1.0)/( total + 1.0)); // add 1.0 to prevent things ever being zero
    }

    return 0;

}


int
write_output(vector<Gene> &genes){
 
    //
    // write a simple 2-column tsv
    // where gene ID -> pvalue
    //

    string outname = "Gene_pvalues.tsv";
    ofstream fh(outname);

    if (!fh.is_open()){
        cerr << "Error: Unable to write " << outname << " in this directory\n";
        exit(1);
    }

    fh << "#Gene\tpvalue\n";
    fh << setprecision(6);

    for (uint i = 0; i < genes.size(); i++){
        fh << genes[i].id << '\t' << genes[i].pval << '\n';
    }

    fh.close();

    return 0;
}

uint
create_uint(const char *arg){

    // if possible, create a number out of this
    // char array

    if (arg == nullptr || *arg == '\0'){
        cerr << "Empty input provided for either --bootstraps, --threads, or --seed\n";
        exit(1);
    }

    const char *end = arg + std::strlen(arg);
    uint val        = 0;
    auto result     = std::from_chars(arg, end, val, 10);

    if (result.ec != std::errc() || result.ptr != end){
        cerr << arg << " is not a valid --bootstraps value\n";
        exit(1);
    }

    return val;
}

void
help(){
    cerr << "Usage: ./bootstrap_gene_fsts -f Avg_gene_fsts.tsv --bootstraps INT [default: 10000] --threads INT [default 1] --seed INT [default 1234] --transcripts [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string infile;
    uint bootstraps       = 10000;
    uint threads          = 1;
    uint seed             = 1234;
    bool transcript_level = false;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-f" && i + 1 < argc){
            infile = argv[i + 1];
        }
        else if (arg == "--transcripts"){
            transcript_level = true;
        }
        else if (arg == "--bootstraps" && i + 1 < argc){
            bootstraps = create_uint(argv[i + 1]);
        }
        else if (arg == "--threads" && i + 1 < argc){
            threads = create_uint(argv[i + 1]);
        }
        else if (arg == "--seed" && i + 1 < argc){
            seed = create_uint(argv[i + 1]);
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

    if (bootstraps < 1){
        cerr << "Bootstraps must be at least 1\n";
        exit(1);
    }

    if (threads < 1){
        cerr << "Threads must be at least 1\n";
        exit(1);
    }

    if (seed < 0){
        cerr << "Seed cannot be negative\n";
        exit(1);
    }

    // read in the data
    vector<Gene> genes;
    load_genes(infile, genes, transcript_level);

    // next, count how often we see each fst
    // value and rank it's frequency based on
    // how many times we see fst scores higher than it
    boostrap(genes, bootstraps, threads, seed);

    // write the output
    write_output(genes);
    
    return 0;
}