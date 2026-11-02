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
using std::fstream;
using std::iostream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::stoi;
using std::sort;
using std::isdigit;
using std::uint8_t;
using std::find;
using std::pair;

//
// parse a bam file to find read alignments that can be circularized to represent the
// mitochondrial genome sequences found in a reference genome
//

typedef unsigned int uint;

//
// Current thresholds used to decided whether a read represents a full mito genome
// TODO: make these command line arguments
//
const double MIN_COVERAGE = 0.90; // fraction of the ref the alignments must cover
const double MAX_SPAN     = 1.5;  // aligned bases / ref length ; fraction of read covering ref (includes wrapping)
const uint   MAX_DIST     = 50;   // dist in bp between two alignments from the same read from the ref's perspective
const uint   READ_OVERLAP = 50;   // dist in bp between two alignments from the same read from the read's perspective

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
    string ref_id;
    char   direction = 'F';
    bool   wraps     = false;
    uint   naligns   = 0;
    uint   start_ref = 0;
    uint   cut       = 0; // position where read was rotated
    int    junction  = 0; // read bases where two alignments wrapped (< 0 == overlap)
    uint   trim_head = 0; // unaligned bases drop from start of read
    uint   trim_tail = 0; // unaligned bases drop from end of read
    uint   nm        = 0; // edit distance reported by bam file
    double coverage  = 0; // fraction of the ref covered

};

class Read {
private:
    vector<Align>          _aligns; // hold onto each alignment
    vector<Reconstruction> _reconstructions; // contains only alignments that can represent an entire mt
    uint8_t               _count = 0;

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

    uint _covered_bases(vector<pair<uint, uint>> &spans){
        
        //
        // count the number of references bases coverd by intervals
        // following half open [start, end) where end is not included
        //

        if (spans.empty()){
            return 0;
        }

        sort(spans.begin(), spans.end());

        uint total     = 0;
        uint cur_start = spans[0].first;
        uint cur_end   = spans[0].second;

        for (uint i = 1; i < spans.size(); i++){
            if (spans[i].first > cur_end){
                //
                // jump to the next interval
                //
                total    += cur_end - cur_start;
                cur_start = spans[i].first;
                cur_end   = spans[i].second;
            }
            else if (spans[i].second > cur_end){
                cur_end = spans[i].second; 
            }
        }

        total += cur_end - cur_start;

        return total;
    }

    string _build(const string &ref_id, char direction, vector<Align> &alns, uint &ref_len){

        //
        // attempt to build a mito genome using these sets of alignments
        // returns a string that best guesses the reason why it failed
        //

        if (alns.empty()){
            return "";
        }
        else if (!this->has_seq()){
            // edge case (very unlikely) where all records are hard clipped
            return "no_sequence"; 
        }

        //
        // step 1: set the alignment order and avoid adding
        //         alignments that cover the same range
        //
        vector<uint> order(alns.size()), chain;

        for (uint i = 0; i < order.size(); i++){
            order[i] = i;
        }

        sort(order.begin(), order.end(), [&alns](uint a, uint b){
            uint len_a = alns[a].read_end - alns[a].read_start;
            uint len_b = alns[b].read_end - alns[b].read_start;
            return len_a != len_b ? len_a > len_b : a < b; // same length pick first
            }
        );

        //
        // now create the starting chain
        //
        for (uint i = 0; i < order.size(); i++){
            const Align &cand = alns[order[i]];
            bool keep         = true;

            for (uint j = 0; j < chain.size(); j++){
                const Align &used = alns[chain[j]];
                int end     = cand.read_end < used.read_end ? cand.read_end : used.read_end;
                int start   = cand.read_start > used.read_start ? cand.read_start : used.read_start;
                int overlap = end - start;
                
                if (overlap >= static_cast<int>(READ_OVERLAP)){
                    keep = false;
                    break;
                }
            } // end of j
            
            if (keep){
                chain.push_back(order[i]);
            }
        } // end of i

        //
        // now sort the chain based on where the read starts
        //
        sort(chain.begin(), chain.end(), [&alns](uint a, uint b){
                if (alns[a].read_start != alns[b].read_start){
                    return alns[a].read_start < alns[b].read_start;
                }
                return alns[a].read_end < alns[b].read_end;
            }
        );

        //
        // step 2: walk along the read. Each alignment has to pick up where the last one
        //         stopped on the reference, or once, jump from the end of the reference
        //         to the beginning
        //
        uint wraps = 0, wrap_idx = 0; // idx in chain of the alignment right after wrap

        for (uint i = 1; i < chain.size(); i++){
            const Align &prev = alns[chain[i - 1]];
            const Align &next = alns[chain[i]];
            int          dist = next.ref_start - prev.ref_end;

            if (dist < 0){
                dist *= -1;
            }
    
            if (dist <= MAX_DIST){
                continue; // continous
            }

            if (prev.ref_end + MAX_DIST >= ref_len && next.ref_start <= MAX_DIST){
                wraps++;
                wrap_idx = i;
            }
            else {
                return "discontinuous"; // the dist b/t alignments exceeded the threshold
            }
        }

        if (wraps > 1){
            return "multiple_wraps"; // this read went >=2 around the mitochondria
        }

        //
        // step 3: is enough the reference covered
        //

        vector<pair<uint, uint>> spans;
        uint nm = 0, span = 0;

        for (uint i = 0; i < chain.size(); i++){
            const Align &aln = alns[chain[i]];
            spans.push_back(std::make_pair(aln.ref_start, aln.ref_end));
            nm += aln.nm;
            span += aln.ref_end - aln.ref_start;
        }

        double cov = static_cast<double>(this->_covered_bases(spans)) / ref_len;
        double spn = static_cast<double>(span) / ref_len;

        if (cov < MIN_COVERAGE){
            return "low_coverage";
        }
        else if (spn > MAX_SPAN){
            return "too_long";
        }

        //
        // step 4: put the read into the orientation of the ref
        //

        string oseq = this->seq, oqua = this->quality;

        //
        // the seq is in a different orientation from the
        // one we are trying to build
        //
        if (this->direction != direction){
            oseq = this->_reverse_complement(oseq);
            oqua = string(oqua.rbegin(), oqua.rend());
        }

        uint head = alns[chain.front()].read_start;
        uint tail = 0;

        for (uint i = 0; i < chain.size(); i++){
            tail = tail > alns[chain[i]].read_end ? tail : alns[chain[i]].read_end;
        }

        //
        // edge case if we didn't parse this correctly or the records are weird
        //
        if (oseq.size() != this->size || oqua.size() != oseq.size() || tail > oseq.size()){
            return "sequence_mismatch";
        }

        //
        // step 5: rotate the sequence & quality scores. If wrapping was found, the
        //         sequence starts where the alignment left of
        //
        Reconstruction rec;

        rec.name      = this->name;
        rec.ref_id    = ref_id;
        rec.direction = direction;
        rec.naligns   = chain.size();
        rec.coverage  = cov;
        rec.nm        = nm;
        rec.trim_head = head;
        rec.trim_tail = oseq.size() - tail;
        rec.start_ref = alns[chain.front()].ref_start + 1;

        uint cut = head;
        if (wraps == 1){
            const Align &before = alns[chain[wrap_idx - 1]];
            const Align &after  = alns[chain[wrap_idx]];

            // update
            cut           = after.read_start;
            rec.wraps     = true;
            rec.junction  = static_cast<int>(after.read_start) - static_cast<int>(before.read_end);
            rec.start_ref = after.ref_start + 1;
        }

        rec.cut = cut;
        rec.seq = oseq.substr(cut, tail - cut) + oseq.substr(head, cut - head);
        rec.qua = oqua.substr(cut, tail - cut) + oqua.substr(head, cut - head);

        this->_reconstructions.push_back(rec);

        return ""; // no reason, this worked

    }


public:

    string name, seq, quality;
    uint   size;
    char   direction;
    vector<string> fail_reasons; // why this read couldn't form a complete mt genome


    Read(string &n){
        this->name = n;
        this->size = 0;
    }

    ~Read(void){
        this->_aligns.clear();
        this->_reconstructions.clear();
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
        this->quality   = quality;
        this->direction = direction;
    }

    uint8_t get_count(void){
        this->_count++;
        return this->_count;
    }

    void reconstruct(unordered_map<string, uint> &lengths){
        //
        // step 1: bin each alignment to their reference alignment
        //
        unordered_map<string, _organizer> _mappings;

        for (uint i = 0; i < this->_aligns.size(); i++){
            Align &aln = this->_aligns[i];

            if (aln.direction == 'F'){
                _mappings[aln.ref_id].fwd_alignments.push_back(aln);
            }
            else {
                _mappings[aln.ref_id].rev_alignments.push_back(aln);
            }
        }

        //
        // step 2: iterate over all possible cases and
        //         attempt to build a mito genome out of it
        //

        string reason;
        uint  ref_len;
        for (auto itr = _mappings.begin(); itr != _mappings.end(); itr++){
            ref_len = lengths[itr->first];
            
            // forward direction
            reason  = this->_build(itr->first, 'F', itr->second.fwd_alignments, ref_len);
            
            // this failed for some reason
            if (!reason.empty()){
                this->fail_reasons.push_back(reason);
            }

            reason = this->_build(itr->first, 'R', itr->second.rev_alignments, ref_len);

            if (!reason.empty()){
                this->fail_reasons.push_back(reason);
            }
        } // end of itr

        //
        // step 3: best explanation first
        //
        sort(this->_reconstructions.begin(), this->_reconstructions.end(), 
            [](const Reconstruction &a, const Reconstruction &b){
                // least mismatches with ref
                if (a.nm != b.nm){
                    return a.nm < b.nm;
                }
                // which rec covered the mtgenome the most
                else if (a.coverage != b.coverage){
                    return a.coverage > b.coverage;
                }
                // sort by alphabetical order
                else if (a.ref_id != b.ref_id){
                    return a.ref_id < b.ref_id;
                }
                // return F first if all things are equal
                return a.direction < b.direction;
            }
        );

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
verify(const string &bamfile, const vector<string> &targets, unordered_map<string, uint> &lengths){
    //
    // helper function to determine whether the targets
    // were part of this alignment file to begin with
    //

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

                if (find(targets.begin(), targets.end(), ref) != targets.end()){
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
            ref = parts[2];
            if (find(targets.begin(), targets.end(), ref) == targets.end()){
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
            aln.id      = query.get_count();

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

    return 0;
}

int 
write_mappings(vector<Read> &mappings, unordered_map<string, uint> &lengths){

    //
    // create the output stream to construct the file
    // containing all the possible full reconstructions
    // in both a tsv file that contains details of the constructions
    // along with a constructed fastq file
    //

    ofstream ofh("Reconstruct_Sequences.tsv");
    ofstream fqfh("Reconstruct_Sequences.fastq");

    uint failed = 0, written = 0;
    unordered_map<string, uint> reasons;

    ofh << "#Read\tRef\tStrand\tRead_length\tN_Alignments\tCoverage\tWraps\t"
        << "Start_Ref\tCut\tJunction\tTrim_Head\tTrim_Tail\tLength\tnm\n";


    for (int i = 0; i < mappings.size(); i++){
        //
        // try to reconstruct the scaffold based on its
        // alignments
        //
        mappings[i].reconstruct(lengths);

        if (!mappings[i].has_reconstructions()){
            failed++;
            // increment the number of reasons
            for (uint j = 0; j < mappings[i].fail_reasons.size(); j++){
                reasons[mappings[i].fail_reasons[j]]++;
            }
            continue; 
        }

        vector<Reconstruction> &reconstructions = mappings[i].get_reconstructions();
        uint                           read_len = mappings[i].size;

        for (uint j = 0; j < reconstructions.size(); j++){
            Reconstruction &rec = reconstructions[j];
            char         strand = rec.direction == 'F' ? '+' : '-';
            string        wraps = rec.wraps ? "Yes" : "No";

            ofh << rec.name       << '\t' << rec.ref_id    << '\t' << strand        << '\t'
                << read_len       << '\t' << rec.naligns   << '\t' << rec.coverage  << '\t'
                << wraps          << '\t' << rec.start_ref << '\t' << rec.cut       << '\t'
                << rec.junction   << '\t' << rec.trim_head << '\t' << rec.trim_tail << '\t'
                << rec.seq.size() << '\t' << rec.nm        << '\n';
        } // end of j

        Reconstruction &best = reconstructions.front();

        fqfh << '@' << best.name << " ref=" << best.ref_id << " start=" 
             << best.start_ref   << '\n'   << best.seq     << "\n+\n" 
             << best.qua         << '\n';
        written++;
    }

    ofh.close();
    fqfh.close();

    cerr << "Wrote " << written  << " reads\n";
    cerr << failed << " out of " << mappings.size() << " reads failed to reconstruct to a mitogenome\n";

    for (auto itr = reasons.begin(); itr != reasons.end(); itr++){
        cerr << '\t' << itr->first << ": " << itr->second << '\n';
    }

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
    unordered_map<string, uint> lengths; // lengths of target
    parse_comma_delim(target_str, targets);

    if (targets.empty()){
        cerr << "Error: Unable to parse " << target_str << " for reference sequence IDs\n";
        exit(EXIT_FAILURE);
    }
    else {
        //
        // verify that the target IDs can be found in the bam file
        //

        verify(bamfile, targets, lengths);
    }

    //
    // step 2: collect all the reads that map to the targets
    //
    vector<Read> mappings;
    collect_alignments(bamfile, mappings, targets);

    //
    // step 3/4: right before writing the collected records
    //           reconstruct a sequence based on its mtgenome
    //           so that it best resembles the target mtgenomes
    //
    write_mappings(mappings, lengths);

    return 0;
}