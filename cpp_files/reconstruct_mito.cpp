#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cctype>
#include <cstdint>

using std::string;
using std::unordered_map;
using std::vector;
using std::fstream;
using std::iostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::stoi;
using std::isdigit;
using std::uint8_t;
using std::find;

//
// parse a bam file to find read alignments that can be circularized to represent the
// mitochondrial genome sequences found in a reference genome
//

typedef unsigned int uint;

struct Align {
    string ref_id;
    string cigar;
    string seq;
    string qual;
    uint   pos;
    char   direction;
    uint8_t id; // id for the alingment, not the read
};

struct Part {
    uint8_t id;
    uint  start; // relative to the reference
    uint  length;

    //
    // we need to define a constructor
    //
    Part(uint8_t id_, uint start_, uint length_) :
        id(id_), start(start_), length(length_) {}
};

struct Reconstruction {
    string name;
    string seq;
    string qua;

};


class Read {

private:
    vector<Align>          _aligns; // hold onto each alignment
    vector<Reconstruction> _reconstructions; // contains only alignments that can represent an entire mt
    uint8_t               _count = 0;


public:

    string name;
    uint   size;

    Read(string &n){
        this->name = n;
    }

    ~Read(void){
        this->_aligns.clear();
        this->_reconstructions.clear();
    }

    void add_align(Align &a){
        this->_aligns.emplace_back(a);
    }

    uint8_t get_count(void){
        this->_count++;
        return this->_count;
    }

    void reconstruct(void){
        //
        // step 1: bin each alignment to their reference alignment
        //
        vector<vector<Align>> _mappings;

        for (uint i = 0; i < this->_aligns.size(); i++){
            Align &aln = this->_aligns[i];
            bool added = false;
            uint j     = 0 ;
            for (; j < _mappings.size(); j++){
                if (_mappings[j].front().ref_id == aln.ref_id){
                    _mappings[j].push_back(aln);
                    added = true;
                    break;
                }
            }
            if (!added){
                _mappings.push_back(vector<Align>());
                _mappings[j].push_back(aln);
            }
        }

        //
        // simple case: all alignments to a single reference
        //

        //
        // the mitochondria may have undergone recombination
        //

    }

    bool has_reconstructions(void){
        return !this->_reconstructions.empty();
    }

    vector<Reconstruction> & get_reconstructions(void){
        return this->_reconstructions;
    }


};

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

int
parse_comma_delim(string &target_str, vector<string> &parts){

    //
    // used once for parsing the targets
    //

    int start = 0, end = 0, n = target_str.size();

    while (end < n){
        if (target_str[end] == ','){
            parts.emplace_back(target_str.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }

    if (start < n){
        parts.emplace_back(target_str.substr(start));
    }

    return 0;
}

int
parse_tabular(string &line, vector<string> &parts){

    // always reset before filling
    parts.clear();

    int start = 0, end = 0, n = line.size();

    while (end < n){
        if (line[end] == '\t'){
            parts.emplace_back(line.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }

    if (start < n){
        parts.emplace_back(line.substr(start));
    }

    return 0;
}

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

bool
check_for_index(const string& bamfile){
    //
    // check if the bam file has an index
    // file in the same path
    //
    return file_exists(bamfile + ".bai") || file_exists(bamfile + ".csi");
}

int
get_sequence_length(const string &cigar){

    //
    // parse the cigar to determine the sequence length
    //

    uint size = 0, p = 0, n = cigar.size(), last_p = 0;
    string cur_str;

    while (p < n){
        char c = cigar[p];
        if (isdigit(c)){
            p++;
        }
        else{
            if (c != 'D'){
                cur_str = cigar.substr(last_p, p - last_p);
                size   += stoi(cur_str);
            }
            last_p = p + 1; // skip this char
            p++;
        }
    }

    return size;
}

bool
is_unmapped(const string &flag){
    //
    // check if this is an unmapped record
    //

    return stoi(flag) & 0x4;
}

bool
is_primary(const string &flag){
    //
    // check if the alignment is primary
    //

    int enc = stoi(flag); // encoding

    if (enc & 0x100 || enc & 0x800){
        return false; }

    return true;
}


bool
reverse_strand(const string &flag){
    //
    // determine if this alignment is on the reverse strand
    //

    return (stoi(flag) & 0x10) != 0;
}

int
verify(const string &bamfile, const vector<string> &targets){
    //
    // helper function to determine whether the targets
    // were part of this alignment file to begin with
    //

    uint found = 0;
    string cmd = "samtools view -H " + bamfile;
    FILE *fh   = popen(cmd.c_str(), "r");
    vector<string> parts;
    string line, region, ref;
    size_t p;
    int c;

    //
    // We will have to build the line
    //
    while ((c = fgetc(fh)) != EOF){
        if (c == '\n'){
            parse_tabular(line, parts);

            if (parts.empty()){
                continue;
            }
            else if (parts.front()[0] != '@'){
                region = parts[1];
                p      = region.find_first_of(':');
                ref    = region.substr(p + 1);
                if (find(targets.begin(), targets.end(), ref) != targets.end()){
                    found++;
                    if (found == targets.size()){
                        break; // done early
                    }
                }
            }
            line.clear();
        }
        else {
            line += static_cast<char>(c);
        }
    }
    pclose(fh);

    if (found != targets.size()){
        cerr << "Error: Not all the target IDs were found in the bam file\n";
        exit(EXIT_FAILURE);
    }

    return 0;

}

int
collect_alignments(const string &bamfile, vector<Read> &mappings, const vector<string> &targets){

    
    //
    // collect all the sequence alignments on the target sequences
    //
    string line, ref;
    vector<string> parts;

    //
    // create a file stream
    //
    string cmd = "samtools view " + bamfile;
    FILE *fh   = popen(cmd.c_str(), "r");

    int c;
    //
    // We will have to build the line
    //
    while ((c = fgetc(fh)) != EOF){
        if (c == '\n'){
            parse_tabular(line, parts);

            if (parts.empty()){
                line.clear();
                continue;
            }
            ref = parts.front();
            if (find(targets.begin(), targets.end(), ref) == targets.end()){
                continue; // not a target chrom
            }

            Align aln;
            string name   = parts[0];
            aln.cigar     = parts[5];
            aln.ref_id    = parts[2];
            aln.seq       = parts[9];
            aln.qual      = parts[10];
            aln.direction = reverse_strand(parts[1]) ? 'R' : 'F';
            aln.pos       = static_cast<uint>(stoi(parts[3]));
            for (int i = 0; i < mappings.size(); i++){
                if (mappings[i].name == name){
                    aln.id = mappings[i].get_count();
                    mappings[i].add_align(aln);
                    name.clear();

                    // check if we need to get the sequence length
                    if (is_primary(parts[1])){
                        mappings[i].size = get_sequence_length(aln.cigar);
                    }
                    break;
                }
            }
            // did we not add this?
            if (!name.empty()){
                Read query(name);
                aln.id = query.get_count();
                query.add_align(aln);
                if (is_primary(parts[1])){
                    query.size = get_sequence_length(aln.cigar);
                }
                mappings.push_back(query);
            }
        parts.clear();
        line.clear();
    }
    else {
        line += static_cast<char>(c);
    }
    }
    pclose(fh);

    if (mappings.empty()){
        cerr << "No alignments were found for the target sequences\n";
        exit(EXIT_FAILURE);
    }

    return 0;
}



int 
write_mappings(vector<Read> &mappings){

    //
    // TODO: This function is incomplete
    // fill finish once reconstruction logic is finished
    //

    //
    // create the output stream to construct the file
    // containing all the possible full reconstructions
    //
    ofstream ofh("Reconstruct_Sequences.tsv");
    ofstream fqfh("Reconstruct_Sequences.fastq");

    uint failed = 0, algn_count;


    for (int i = 0; i < mappings.size(); i++){
        //
        // try to reconstruct the scaffold based on its
        // alignments
        //
        mappings[i].reconstruct();

        if (!mappings[i].has_reconstructions()){
            failed++;
            continue; // was too fragmented
        }

        vector<Reconstruction> &reconstructions = mappings[i].get_reconstructions();

    }

    ofh.close();
    fqfh.clear();

    cerr << failed << " out of " << mappings.size() << " scaffolds failed to reconstruct\n";

    return 0;
}

void
help(void){
    cerr << "./reconstruct_mito.cpp -b input.bam -t mito1,mito2,etc..\n"
         << "\t-b: Bam file to parse & obtained reads from [required]\n"
         << "\t-t: Comma separated list of mitochondrial references in bam file [required]\n"
         << "\t-h: Help message\n";
    exit(EXIT_SUCCESS);
}

int
main(int argc, char *args[]){

    string bamfile, arg, target_str;

    for (int i = 1; i < argc; i++){
        arg = args[i];
        if (arg == "-b" && i + 1 < argc){
            bamfile = args[i + 1];
        }
        else if (arg == "-t" && i + 1 < argc){
            target_str = args[i + 1];
        }
        else if (arg == "-h"){
            help();
        }

    }

    if (bamfile.empty()){
        cerr << "No bam file was provided. Program requires a bam file (-b bam)\n";
        exit(EXIT_FAILURE);
    }
    else if (file_exists(bamfile) == false){
        cerr << "Could not locate " << bamfile << '\n';
        exit(EXIT_FAILURE);
    }
    else if (check_for_index(bamfile) == false){
        cerr << "Could not locate a .bai or .csi index file for " << bamfile << '\n';
        exit(EXIT_FAILURE);
    }

    if (target_str.empty()){
        cerr << "-t requires an input string\n";
        exit(EXIT_FAILURE);
    }

    //
    // step 1: collect the ids of the target reference mito genomes
    //
    vector<string> targets;
    parse_comma_delim(target_str, targets);
    if (targets.empty()){
        cerr << "Error: Unable to parse " << target_str << " for reference sequence IDs\n";
        exit(EXIT_FAILURE);
    }
    else {
        //
        // verify that the target IDs can be found in the bam file
        //
        verify(bamfile, targets);
    }

    //
    // step 2: collect all the reads that map to the targets
    //
    vector<Read> mappings;
    collect_alignments(bamfile, mappings, targets);

    //
    // step 2: write these out to a file
    //
    write_mappings(mappings);

    return 0;
}