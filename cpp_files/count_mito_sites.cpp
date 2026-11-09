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
using std::pair;
using std::cout;
using std::cerr;
using std::stoi;
using std::sort;
using std::lower_bound; // binary search
using std::unique;
using std::toupper;
using std::isdigit;
using std::uint8_t;
using std::getline;

//
// given lists of fixed differences [2 column tsv files] between two mitotypes, count the number of
// positions that match one of the provided mitotypes
//

typedef unsigned int uint; // for depth calculations
const char nucleotides[5] = {'A', 'C', 'G', 'T', 'N'};

struct Site {
    uint pos;
    char nuc;

    Site(uint pos_, char nuc_) :
        pos(pos_), nuc(nuc_) {}
};

struct SiteDepth {
    uint counts[5] = {0, 0, 0, 0 ,0};
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
};

class Read {
private:
    vector<Align>  _aligns; // hold onto each alignment

    struct _organizer {
        vector<Align> fwd_alignments;
        vector<Align> rev_alignments;
    };

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
        this->seq.clear();
    }

    bool has_seq(void){
        return !this->seq.empty();
    }

    bool has_full_seq(void){
        return !this->seq.empty() && this->seq.size() == this->size;
    }

    void set_seq(const string &seq, char direction){
        
        if (this->has_seq()){
            return; // seq already added
        }

        this->seq       = seq;
        this->direction = direction;
    }

    void add_align(Align &a){
        this->_aligns.emplace_back(a);
    }

    vector<Align> & get_alignments(void){
        return this->_aligns;
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
parse_map(const string &map_file, unordered_map<string, string> &targets, vector<string> &order){

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
        order.push_back(target);
    }

    fh.close();

    if (targets.empty()){
        cerr << "No targets were found in the mapping file\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loaded " << targets.size() << " targets\n";

    return 0;

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
                static_cast<uint>(stoi(parts[0]));
            }
            catch (const std::invalid_argument &error){
                cerr << "Error: Column 1 must be an integer. Found " << parts[0]
                     << " in " << itr->second << '\n';
                exit(EXIT_FAILURE);
            }
            if (parts[1].size() != 1){
                cerr << "Error: Column 2 must be a single nucleotide. Found " << parts[1]
                     << " in " << itr->second << '\n';
                exit(EXIT_FAILURE);
            }
            else if (parts[1] != "A" && parts[1] != "T" && parts[1] != "G" && parts[1] != "C"){
                cerr << "Error: Invalid character in column 2 in sites file. Found " << parts[1]
                     << " in " << itr->second << '\n';
                exit(EXIT_FAILURE);
            }
            
            uint pos = static_cast<uint>(stoi(parts[0]));
            char nuc = parts[1][0];

            sites[target].emplace_back(pos, nuc);
        } // end of while

        fh.close();

        // sort the sites just to be safe
        sort(sites[target].begin(), sites[target].end(), []
            (const Site &a, const Site &b){
                return a.pos < b.pos;}
            );
    }

    cerr << " Done\n";

    for (auto itr = sites.begin(); itr != sites.end(); itr++){
        cerr << '\t' << itr->first << ": " << itr->second.size() << " sites\n";
    }

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
reverse_strand(const string &flag){
    //
    // determine if this alignment is on the reverse strand
    //

    return (stoi(flag) & 0x10) != 0;
}

string
reverse_complement(const string &seq){
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
            if (targets.count(ref) == 0 || is_unmapped(parts[1]) || parts[5] == "*"){
                line.clear();
                continue; // not a target chrom or not aligned
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

            for (uint i = 0; i < mappings.size(); i++){
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
                //
                query.set_seq(parts[9], aln.direction);
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
        cerr << "Found " << mappings.size() << " reads mapped to target mito genomes\n";
    }

    return 0;
}

uint
base_index(char b){
    // helper function to get the idx
    // found in nucleotides
    switch (toupper(b)){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}

int
scan_alignment(const Align &aln, const string &seq, const vector<Site> &markers, vector<uint> &hits, vector<pair<uint, uint>> &bases){
    //
    // helper function to walk through an alignment's cigar and read's
    // sequence and fill a hits array so that overlapping alignments
    // does not contribute more than once. The bases array holds every
    // (site, nuc) from this read, whether it matches 
    // or not
    //

    const string &cigar = aln.cigar;
    uint cur  = 0, prev = 0, n = cigar.size();
    uint rpos = aln.ref_start; // 0-based pos on reference
    uint qpos = 0;             // 0-based pos on read

    while (cur < n){
        char c = cigar[cur];
        if (isdigit(c)){
            cur++;
        }
        else {
            uint len = stoi(cigar.substr(prev, cur - prev));

            if (c == 'M' || c == '=' || c == 'X'){
                // this block is aligned to the referenc and is 1 based.
                // we need to calculate where we are on the reference then
                // see which markers are from here and scan to the end
                // of the align block
                //
                uint block_start = rpos + 1;
                uint block_end   = rpos + len;

                auto itr = lower_bound(markers.begin(), markers.end(), block_start,
                                      [](const Site &s, uint p){ return s.pos < p;});

                // do not go past the markers array or block length
                for (; itr != markers.end() && itr->pos <= block_end; itr++){
                    uint idx = qpos + (itr->pos - 1 - rpos);
                    if (idx < seq.size()){
                        // this read is covered
                        bases.emplace_back(static_cast<uint>(itr - markers.begin()), base_index(seq[idx]));
                        if (seq[idx] == itr->nuc){
                            hits.push_back(itr - markers.begin());
                        }
                    }
                }

                rpos += len;
                qpos += len;
            }
            else if (c == 'I' || c == 'S' || c == 'H'){
                qpos += len; // consumes only the read
            }
            else if (c == 'D' || c == 'N'){
                rpos += len; // consumes only the reference/target
            }
            prev = cur + 1;
            cur++;
        }
    }

    return 0;
}

int
count_sites(vector<Read> &mappings, unordered_map<string, vector<Site>> &sites, unordered_map<string, vector<SiteDepth>> &depth){

    //
    // this function will loop through all the aligned reads and for each read
    // count the number of sites that match the targets in the sites mapping
    // where the key is the target chrom ID. Overlapping alignments are 
    // taken account of. We also calculate depth for each target site
    //

    for (auto itr = sites.begin(); itr != sites.end(); itr++){
        depth[itr->first].assign(itr->second.size(), SiteDepth());
    }

    for (uint i = 0; i < mappings.size(); i++){
        Read &read = mappings[i];

        //
        // each read will be given a default count of 0
        //
        for (auto itr = sites.begin(); itr != sites.end(); itr++){
            read.counts[itr->first] = 0;
        }

        if (!read.has_full_seq()){
            // this read likely came from a hard clipped alignment
            continue;
        }

        const string &seq = read.seq;
        string    rev_seq;
        unordered_map<string, vector<uint>> hits; // target -> sites with matching nuc
        unordered_map<string, vector<pair<uint, uint>>> base_counts; // target -> sites with an aligned base
        vector<Align> &algns = read.get_alignments();

        for (uint j = 0; j < algns.size(); j++){
            const Align &aln = algns[j];
            auto         itr = sites.find(aln.ref_id);

            // won't happen but just good practice
            if (itr == sites.end()){
                continue;
            }

            if (aln.direction == read.direction){
                scan_alignment(aln, seq, itr->second, hits[aln.ref_id], base_counts[aln.ref_id]);
            }
            else {
                if (rev_seq.empty()){
                    // compute once
                    rev_seq = reverse_complement(seq);
                }
                scan_alignment(aln, rev_seq, itr->second, hits[aln.ref_id], base_counts[aln.ref_id]);
            }
        } // end of algns loop

        // deduplicate sites
        for (auto itr = hits.begin(); itr != hits.end(); itr++){
            vector<uint> &h = itr->second;

            // only count each site once
            sort(h.begin(), h.end());
            h.erase(unique(h.begin(), h.end()), h.end());
            read.counts[itr->first] = h.size();
        }

        // repeat for depth
        for (auto itr = base_counts.begin(); itr != base_counts.end(); itr++){
            vector<pair<uint, uint>> &b = itr->second;
            sort(b.begin(), b.end());
            b.erase(unique(b.begin(), b.end()), b.end());
            
            //
            // one read -> one increment
            //
            for (uint k = 0; k < b.size(); k++){
                depth[itr->first][b[k].first].counts[b[k].second]++;
            }
        }

        
    }

    return 0;
}

int
write_counts(vector<Read> &mappings, vector<string> &order){

    //
    // helper function to write the tsv file
    //

    ofstream ofh("Mito_Site_Counts.tsv");

    if (!ofh.is_open()){
        cerr << "Error: Unable to create the counts file\n";
        exit(EXIT_FAILURE);
    }

    ofh << "#Read";

    for (uint i = 0; i < order.size(); i++){
        ofh << '\t' << order[i];
    }
    ofh << '\n';

    // for printing stats
    uint mixed = 0, empty = 0, written = 0, specific = 0, skipped = 0;

    for (uint i = 0; i < mappings.size(); i++){
        
        Read &read = mappings[i];

        if (!read.has_full_seq()){
            skipped++;
            continue;
        }

        uint targets = 0;
        ofh << read.name;

        for (uint j = 0; j < order.size(); j++){
            uint count = read.counts[order[j]];
            ofh << '\t' << count;

            if (count > 0){
                targets++;
            }
        }

        ofh << '\n';
        written++;

        if (targets > 1){
            mixed++;
        }
        else if (targets == 1){
            specific++;
        }
        else {
            empty++;
        }
    }

    ofh.close();

    cerr << "Wrote " << written  << " reads to Mito_Site_Counts.tsv\n"
         << '\t'     << specific << " reads had sites to only one target\n"
         << '\t'     << mixed    << " reads had sites for more than one target\n"
         << '\t'     << empty    << " reads did not have any sites\n"
         << '\t'     << skipped  << " reads were skipped because their full sequence was not in the bam file\n";

    return 0;
}

int
write_depth(unordered_map<string, vector<Site>> &sites, unordered_map<string, vector<SiteDepth>> &depth, vector<string> &order){

    //
    // one row per allele seen at a site: the number of reads with that
    // nucleotide. 0 counts are ignored
    //

    ofstream ofh("Mito_Site_Depths.tsv");

    if (!ofh.is_open()){
        cerr << "Error: Unable to create the depth file\n";
        exit(EXIT_FAILURE);
    }

    ofh << "#Target\tPos\tTarget_Nuc\tAllele\tCount\n";

    for (uint i = 0; i < order.size(); i++){

        const vector<Site>       &tsites = sites[order[i]]; // target sites
        const vector<SiteDepth> &tdepths = depth[order[i]];

        for (uint j = 0; j < tsites.size(); j++){
            for (uint k = 0; k < 5; k++){
                if (tdepths[j].counts[k] == 0){
                    continue;
                }
                ofh << order[i]             << '\t' << tsites[j].pos  << '\t' 
                    << tsites[j].nuc        << '\t' << nucleotides[k] << '\t' 
                    << tdepths[j].counts[k] << '\n';
            }
        }
    }

    ofh.close();

    cerr << "Wrote site depths to Mito_Site_Depths.tsv\n";

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
    vector<string> order;
    parse_map(map_file, targets, order);
    
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
    unordered_map<string, vector<SiteDepth>> depth; // target -> depth of each site
    count_sites(mappings, sites, depth);

    //
    // step 6: write the results
    //
    write_counts(mappings, order);
    write_depth(sites, depth, order);

    return 0;
}