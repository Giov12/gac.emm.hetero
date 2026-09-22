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
using std::stoi;
using std::setprecision;


//
// this code is hard coded for a specific set of samples
// use and modification of this code is free, but as is,
// this code will not be appropriate for your data
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
pop_of_interest(const string &pop){
    //
    // check if this is a population of interest
    //
    if (pop.size() < 8){
        return false; // 8 is the length of the prefix
    }
    return pop.substr(0, 8) == "bear_paw";
}

int
parse_table(const string &table){

    //
    // this function will be the main work horse for this
    // code
    //

    bool gzipped = is_compressed(table);
    gzFile gz_fh = NULL;
    ifstream txt_fh;
    ofstream ofh("Heterozygous_candidates.tsv");

    if (gzipped){
        gz_fh = gzopen(table.c_str(), "rb");
        if (gz_fh == NULL){
            cerr << "Error: could not open " << table << '\n';
            exit(1);
        }
    }
    else {
        txt_fh.open(table);
        if (!txt_fh.is_open()){
            cerr << "Error: could not open " << table << '\n';
            exit(1);  
        }
    }
  
    vector<string> parts;
    vector<uint> ref_pops, het_pops;
    vector<string> samples;
    string   line, pop, column;
    bool     eof = false; // default value

    // some counters
    uint total_sites = 0, line_cnt = 0, total_markers = 0;

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

        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        // this should not be triggerd
        if (line.empty()){
            continue;
        }

        parse_tabular(line, parts);

        if (line_cnt == 1){
            // header line
            for (uint i = 2; i < parts.size(); i++){
                if (pop_of_interest(parts[i])){
                    het_pops.push_back(i);
                }
                else {
                    ref_pops.push_back(i);
                }
            }
            ofh << line << '\n';

            continue;
        }

        total_sites++;

        //
        // now we need to check the reference populations
        // to determine if they are all homozygous
        // first, check that hets are at 0 & if true,
        // ensure that there are at least counts for homo
        // i.e., this site isn't completely missing
        //
        bool candidate = true;
        for (uint i = 0; i < ref_pops.size(); i++){
            Entry e;
            column = parts[ref_pops[i]];
            parse_column(column, e);
            if (e.hetero > 0 || (e.homo_alt == 0 && e.homo_ref == 0)){
                candidate = false;
                break;
            }
        }

        // not useful for our study
        if (!candidate){
            continue;
        }
        
        // now to verify that the other ref populations are heterozygous
        candidate = true;
        for (uint i = 0; i < het_pops.size(); i++){
            Entry e;
            column = parts[het_pops[i]];
            parse_column(column, e);
            if (e.hetero == 0 || e.homo_alt > 0 || e.homo_ref > 0){
                candidate = false;
                break;
            }
        }

        // this marker survived all tests
        if (candidate){
            total_markers++;

            // now to keep this line
            ofh << line << '\n';

        }
    } // end of file parsing

    if (gzipped){
        gzclose(gz_fh);
    }
    else {
        txt_fh.close();
    }

    ofh.close();

    if (total_sites == 0){
        cerr << "Found no records in this file\n";
        return 1;
    }

    double p = ((double)total_markers/ (double)total_sites) * 100.0;
    cerr << std::fixed << setprecision(2);
    cerr << "Scanned " << total_sites << " sites and found " <<  total_markers << " (" << p << "%) candidates\n";

    return 0;
}

void
help(){
    cerr << "Usage: ./find_hwe_sites.cpp -t merged_hwe.tsv.gz\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string table;
    
    // expect at least a single argument
    if (argc < 2){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-t" && i + 1 < argc){
            table = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (table.empty()){
        help();
    }

    if (!file_exists(table)){
        cerr << "Unable to find " << table << '\n';
        exit(1);
    }

    // parse the table
    parse_table(table);

    return 0;
}