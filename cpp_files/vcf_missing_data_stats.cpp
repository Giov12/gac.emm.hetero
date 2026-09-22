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
// of missing data for each sample in a vcf file
//

typedef unsigned int uint;

struct Sample {
    string sample;
    uint   count; // missing count
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
parse_vcf(const string &vcf){

    //
    // this function will be the main work horse for this
    // code
    //

    bool gzipped = is_compressed(vcf);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    if (gzipped){
        gz_fh = gzopen(vcf.c_str(), "rb");
        if (gz_fh == NULL){
            cerr << "Error: could not open " << vcf << '\n';
            exit(1);
        }
    }
    else {
        txt_fh.open(vcf);
        if (!txt_fh.is_open()){
            cerr << "Error: could not open " << vcf << '\n';
            exit(1);  
        }
    }
    //
    // sample -> number of missing sites for this sample
    //
    

    vector<string> parts;
    vector<Sample> samples;
    string   line, sample;
    bool     eof = false; // default value

    //
    // calculate heterozygosity at every exonic
    // site for each population
    //
    uint total_sites = 0, total_with_all = 0;

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
                    samples.push_back({sample, 0});
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
        total_sites++;

        uint missing = 0;
        for (uint i = 9; i < parts.size(); i++){
            uint sample_idx = i - 9;
            string geno     = parts[i];
            geno            = geno.substr(0, geno.find(':')); 
            if (geno.find('.') != string::npos){ // missing data
                samples[sample_idx].count++;
                missing++;
            }
        }
        if (missing == 0){
            total_with_all++;
        }

    } // end of file parsing

    if (gzipped){
        gzclose(gz_fh);
    }
    else {
        txt_fh.close();
    }

    if (total_sites == 0){
        cerr << "No variants sites found in " << vcf << '\n';
        return 1;
    }
    double all_geno = ((double)total_with_all / (double)total_sites) * 100.0;

    cerr << std::fixed << std::setprecision(2);
    cerr << "Total number of sites " << total_sites << '\n'
         << "Total number of sites with all genotyped " << total_with_all << " (" << all_geno << "%)\n";
    
    for (uint i = 0; i < samples.size(); i++){
        double miss = ((double)samples[i].count / (double)total_sites) * 100.0;
        cerr << "Sample " << samples[i].sample << " is missing " << miss << "% of sites\n";
    }

    return 0;
}

void
help(){
    cerr << "Usage: ./vcf_missing_data_stats -v vcf_file\n";
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