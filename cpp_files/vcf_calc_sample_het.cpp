#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <zlib.h>
#include <iomanip>

using std::string;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
using std::cout;
using std::setprecision;


//
// some code I wrote just to get an idea on the amount
// heterozygous sites per sample
//

typedef unsigned int uint;

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

int
parse_vcf(const string &vcf){

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
    string   line, sample;
    bool     eof = false; // default value

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
    cerr << "Usage: ./vcf_calc_sample_het -v vcf_file\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf;
    
    // expect at least a single argument
    if (argc < 2){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = string(argv[i + 1]);
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

    // parse the vcf file
    parse_vcf(vcf);

    return 0;
}