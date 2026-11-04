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
#include <utility>

using std::string;
using std::unordered_map;
using std::vector;
using std::ifstream;
using std::iostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::stoi;
using std::sort;
using std::isdigit;
using std::uint8_t;
using std::getline;

//
// given lists of fixed differences [2 column tsv files] between two mitotypes, count the number of
// positions that match one of the provided mitotypes
//

typedef unsigned int uint;

struct Site {
    int16_t pos;
    char    nuc;
};
struct Align {
    string ref_id;
    string cigar;
    uint   pos;
    uint   ref_start;
    uint   ref_end;
    uint   read_start;
    uint   read_end;
    char   direction;
    bool   hard_clipped;
    uint    nm; // edit distance from bam file
    uint8_t id; // id for the alingment, not the read  
};


class Read {
private:
    vector<Align>  _aligns; // hold onto each alignment

    struct _organizer {
        vector<Align> fwd_alignments;
        vector<Align> rev_alignments;
    };

    string _reverse_complement(const string &seq){
        //
        // return the reverse compliment of this sequence
        //

        const uint n = seq.size();
        string rc(n, 'N'); // string of size n full of N's

        for (uint i = 0; i < n; i++){
            switch (seq[n - i - 1]) {
                case 'A':
                    rc[i] = 'T';
                    break;
                case 'C':
                    rc[i] = 'G';
                    break;
                case 'G':
                    rc[i] = 'C';
                    break;
                case 'T':
                    rc[i] = 'A';
                    break;
                case 'a':
                    rc[i] = 'T';
                    break;
                case 'c':
                    rc[i] = 'G';
                    break;
                case 'g':
                    rc[i] = 'C';
                    break;
                case 't':
                    rc[i] = 'A';
                    break;
                default:
                    rc[i] = 'N'; // unknown char
                    break;
            }
        }

        return rc;
    }

public:

    string                   name, seq;
    uint                     size;
    char                     direction;
    unordered_map<string, uint> counts;


    Read(string &n){
        this->name = n;
        this->size = 0;
    }

    ~Read(void){
        this->_aligns.clear();
    }

    bool has_seq(void){
        return !this->seq.empty();
    }

    void add_align(Align &a){
        this->_aligns.emplace_back(a);
    }

    void set_seq(const string &seq, const string &quality, char direction){

        if (this->has_seq()){
            return; // sequence already added
        }

        this->seq       = seq;
        this->direction = direction;
    }

};

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
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

int
parse_cigar(Align &aln){
    
    //
    // walk through each operation & get the coordinates
    // both on the reference and the read
    //

    const string &cigar  = aln.cigar;
    uint cur   = 0, prev = 0,     n = cigar.size();
    uint lead  = 0, qlen = 0, rlen = 0;
    bool algnd = false; // did we pass the leading clip

    aln.hard_clipped = false; // until proven

    while (cur < n){
        char c = cigar[cur];
        if (isdigit(c)){
            cur++;
        }
        else {
            uint len = stoi(cigar.substr(prev, cur - prev));
            if (c == 'S' || c == 'H'){
                if (!algnd){
                    // leading clip, trailing clip is implied by the read length
                    lead += len; 
                }
                if (c == 'H'){
                    aln.hard_clipped = true;
                }
            }
            else {
                algnd = true;
                //
                // see https://samtools.github.io/hts-specs/SAMv1.pdf
                //
                if (c == 'M' || c == 'I' || c == '=' || c == 'X'){
                    qlen += len; // consumes the read
                }
                if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X'){
                    rlen += len; // consumes the reference
                }
            }
            prev = cur + 1;
            cur++;
        }
    }

    aln.read_start = lead;
    aln.read_end   = lead + qlen;
    aln.ref_start  = aln.pos - 1;
    aln.ref_end    = aln.ref_start + rlen;

    return 0;
}

int
parse_map(const string &map_file, unordered_map<string, string> &targets){

    //
    // read in the tsv file to ensure that the target files
    // exist
    //

    ifstream fh(map_file);

    if (!fh.is_open()){
        cerr << "Error: Unable to open " << map_file << '\n';
        exit(EXIT_FAILURE);
    }

    string line, target, path;
    vector<string> parts;

    while (getline(fh, line)){
        parse_tabular(line, parts);

        if (parts.size() != 2){
            cerr << "Error: Invalid format in mapping file. Expected 2 columns."
                 << " Offending line: " << line << '\n';
            exit(EXIT_FAILURE);
        }
        target = parts[0];
        path   = parts[1];

        if (file_exists(path) == false){
            cerr << "Error: Could not locate " << path << '\n';
            exit(EXIT_FAILURE);
        }
        targets[target] = path;
    }

    fh.close();

    if (targets.empty()){
        cerr << "No targets were found in the mapping file\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loaded " << targets.size() << " targets\n";

};

int
load_sites(unordered_map<string, string> &targets, unordered_map<string, vector<Site>> &sites){

    //
    // helper function to load all the fixed sites for each of the
    // mitotypes
    //

    cerr << "Loading sites to each target...";

    for (auto itr = targets.begin(); itr != targets.end(); itr++){

        string target = itr->first;
        ifstream fh(itr->second);

        if (!fh.is_open()){
            cerr << "Error: Unable to open " << itr->second << '\n';
            exit(EXIT_FAILURE);
        }
        string line;
        vector<string> parts;
        while (getline(fh, line)){
            parse_tabular(line, parts);
            
            if (parts.size() != 2){
                cerr << "Error: Invalid format in " << itr->second << '\n'
                     << "Expected 2 columns. Offending line:\n" << line << '\n';
                exit(EXIT_FAILURE);
            }

            //
            // test for valid inputs
            //
            try {
                static_cast<int16_t>(stoi(parts[0]));
            }
            catch (const std::invalid_argument){
                cerr << "Error: Column 1 must be an integer. Found " << parts[0]
                     << " in " << itr->second << '\n';
                exit(EXIT_FAILURE);
            }
            if (parts[1].size() != 1){
                cerr << "Error: Column 2 must be a single nucleotide. Found " << parts[1]
                     << " in " << itr->second << '\n';
                exit(EXIT_FAILURE);
            }
            else if (parts[1] != "A" || parts[1] != "T" || parts[1] != "G" || parts[1] != "C"){
                cerr << "Error: Invalid character in column 2 in sites file. Found " << parts[1]
                     << " in " << itr->second << '\n';
                exit(EXIT_FAILURE);
            }
            
            int16_t pos = static_cast<int16_t>(stoi(parts[0]));
            char    nuc = parts[1][0];

            sites[target].emplace_back(pos, nuc);
        }
    }

    cerr << " Done\n";

    return 0;
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
            if (c == 'M' || c == 'I' || c == 'S' || c == 'H' || c == '=' || c == 'X'){
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
verify(const string &bamfile, const unordered_map<string, string> &targets, unordered_map<string, uint> &lengths){
    //
    // helper function to determine whether the targets
    // were part of this alignment file to begin with
    //

    cerr << "Verifying " << bamfile << '\n';

    uint found = 0;
    string cmd = "samtools view -H " + bamfile;
    FILE *fh   = popen(cmd.c_str(), "r");
    vector<string> parts;
    string line, region, ref;
    uint ref_len;
    int c;

    //
    // We will have to build the line
    //
    while ((c = fgetc(fh)) != EOF){
        if (c == '\n'){
            parse_tabular(line, parts);

            if (parts.size() > 1 && parts[0] == "@SQ"){
                ref.clear();
                ref_len = 0;

                for (uint i = 1; i < parts.size(); i++){
                    if (parts[i].compare(0, 3, "SN:") == 0){
                        ref = parts[i].substr(3);
                    }
                    else if (parts[i].compare(0, 3, "LN:") == 0){
                        ref_len = static_cast<uint>(stoi(parts[i].substr(3)));
                    }
                }

                if (targets.count(ref) > 0){
                    found++;
                    lengths[ref] = ref_len;
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
    else {
        cerr << "Found all target mt genomes in " << bamfile << '\n';
    }

    return 0;
}

int
collect_alignments(const string &bamfile, vector<Read> &mappings, const unordered_map<string, string> &targets){
    
    //
    // collect all the sequence alignments on the target sequences
    //
    string line, ref;
    vector<string> parts;

    cerr << "Parsing " << bamfile << '\n';

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
            ref = parts[2];
            if (targets.count(ref) == 0){
                line.clear();
                continue; // not a target chrom
            }

            Align aln;
            string name   = parts[0];
            aln.cigar     = parts[5];
            aln.ref_id    = parts[2];
            aln.direction = reverse_strand(parts[1]) ? 'R' : 'F';
            aln.pos       = static_cast<uint>(stoi(parts[3]));
            uint idx      = mappings.size();

            //
            // parse the cigar & add the other attributes (e.g., ref_start)
            //
            parse_cigar(aln);

            //
            // check if the nm flag was set
            //
            aln.nm = 0;
            for (uint i = 11; i < parts.size(); i++){
                // equivalent to str.startswith()
                if (parts[i].compare(0, 5, "NM:i:") == 0){
                    aln.nm = static_cast<uint>(stoi(parts[i].substr(5)));
                    break;
                }
            }

            for (int i = 0; i < mappings.size(); i++){
                if (mappings[i].name == name){
                    idx = i;
                    break;
                }
            }
            // it was not added
            if (idx == mappings.size()){
                mappings.emplace_back(name);
            }
            Read &query = mappings[idx];

             // set once
            if (query.size == 0){ 
                query.size = get_sequence_length(aln.cigar);  
            }

            if (!aln.hard_clipped && parts[9] != "*"){
                //
                // parts[9]  == query sequence
                // parts[10] == query sequence quality scores
                //
                query.set_seq(parts[9], parts[10], aln.direction);
            }
            query.add_align(aln);
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
    else {
        cerr << "Found " << mappings.size() << " alignments mapped to target mito genomes\n";
    }

    return 0;
}

void
help(void){
    cerr << "./count_mito_sites -b input.bam -m map_file.tsv\n"
         << "\t-b: Bam file to parse & obtained reads from [required]\n"
         << "\t-m: A two column tsv that lists the target chrom and path to input file [required]\n"
         << "\t-h: Help message\n";
    exit(EXIT_SUCCESS);
}

int
main(int argc, char *args[]){

    string bamfile, arg, map_file;

    for (int i = 1; i < argc; i++){
        arg = args[i];
        if (arg == "-b" && i + 1 < argc){
            bamfile = args[i + 1];
        }
        else if (arg == "-m" && i + 1 < argc){
            map_file = args[i + 1];
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

    if (map_file.empty()){
        cerr << "-m is required\n";
        exit(EXIT_FAILURE);
    }
    else if (file_exists(map_file) == false){
        cerr << "Could not locate " << map_file << '\n';
        exit(EXIT_FAILURE);
    }

    //
    // step 1: collect the ids of the targets and the file paths
    //
    unordered_map<string, string> targets;
    unordered_map<string, vector<Site>> sites;
    unordered_map<string, uint> lengths; // lengths of target
    parse_map(map_file, targets);
    
    //
    // step 2: verify that these target chroms are in the bam file
    //
    verify(bamfile, targets, lengths);

    //
    // step 3: collect the sites
    //
    load_sites(targets, sites);
    
    //
    // step 4: collect all the reads that map to the targets
    //
    vector<Read> mappings;
    collect_alignments(bamfile, mappings, targets);

    //
    // step 5: now process each read
    //

    return 0;
}