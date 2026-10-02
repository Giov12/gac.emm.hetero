#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include <iomanip>

using std::string;
using std::fstream;
using std::ofstream;
using std::ifstream;
using std::unordered_map;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::stod;
using std::sort;
using std::setprecision;

typedef unsigned int uint;

struct Exon {
    uint start;
    uint end;
};

struct Site {
    uint   pos;
    double value;
};

struct Attribute {
    string key;
    string val;
};

struct Transcript {
    string id;
    double avg;
    uint length = 0;
    vector<Exon> exons;
    vector<double> scores;
};

class Gene {

private:
    // key -> transcript ID
    unordered_map<string, Transcript> _transcripts;
    bool _coding    = false;
    bool _has_score = false;

public:
    string id;
    string name;
    uint   start;
    uint   end;

    Gene (string id_, string name, uint start, uint end){
        this->id    = id_;
        this->name  = name;
        this->start = start;
        this->end   = end;
    };

    ~Gene(){
        this->_transcripts.clear();
    }

    void add_exon(Exon exon, const string &transcript_id){
        //
        // just add the exon for now
        // & then we will sort them later
        //
        this->_transcripts[transcript_id].length += (exon.end - exon.start + 1);
        this->_transcripts[transcript_id].exons.push_back(exon);
    }

    void add_site(Site &site){
        //
        // add the score to each transcript it overlaps
        //
        
        // should not happen
        if (this->_transcripts.empty()){
            return;
        }

        this->_has_score = true; // overwriting is fine

        for (auto itr = this->_transcripts.begin(); itr != this->_transcripts.end(); itr++){
            vector<Exon> &exons = itr->second.exons;
            for (uint i = 0; i < exons.size(); i++){
                if (exons[i].start <= site.pos && site.pos <= exons[i].end){
                    itr->second.scores.push_back(site.value);
                }
            }
        }
    }

    void sort_exons(void){
        //
        // for each set of exons,
        // sort by basepair
        //

        // edge-case
        if (this->_transcripts.empty()){
            return;
        }

        for (auto itr = this->_transcripts.begin(); itr != this->_transcripts.end(); itr++){
            vector<Exon> &exons = itr->second.exons;

            sort(exons.begin(), exons.end(), []
                (const Exon &exon1, const Exon &exon2){
                    if (exon1.start == exon2.start){
                        return exon1.end < exon2.end;
                    }
                    return exon1.start < exon2.start;
                }
            );
        }
    }

    void merge_exons(void){
        
        //
        // create a single Transcript that will
        // encompass all of this gene's exons
        //

        // nothing to do hear
        if (this->_transcripts.empty()){
            return;
        }

        vector<Exon> _all_exons;
        string tname;

        for (auto itr = this->_transcripts.begin(); itr != this->_transcripts.end(); itr++){
            tname               = itr->first; // keep just in case
            vector<Exon> &exons = itr->second.exons;
            for (uint i = 0; i < exons.size(); i++){
                _all_exons.push_back(exons[i]);
            }
            exons.clear();
        }

        //
        // sort to then just go exon by exon
        //
        sort(_all_exons.begin(), _all_exons.end(), []
            (const Exon &exon1, const Exon &exon2){
                if (exon1.start == exon2.start){
                    return exon1.end < exon2.end;
                }
                return exon1.start < exon2.start;
            }
        );

        vector<Exon> resolved;
        resolved.reserve(_all_exons.size());

        resolved.push_back(_all_exons.front());
        int i = 1, count = _all_exons.size();

        while (i < count){
            Exon &prev = resolved.back();
            Exon &next = _all_exons[i];

            // is there overlap?
            if (next.start <= prev.end){
                if (next.end > prev.end){ // merge if true
                    prev.end = next.end;
                }
            }
            else {
                resolved.push_back(next);
            }
            i++;
        }
        
        // clear up transcript map
        string name = this->_transcripts.size() == 1 ? tname : "merged_transcripts";
        Transcript t;
        t.id    = name;
        t.exons = resolved;
        
        for (uint i = 0; i < resolved.size(); i++){
            t.length += resolved[i].end - resolved[i].start + 1;
        }

        this->_transcripts.clear();
        this->_transcripts[name] = t;

    }

    vector<Transcript*> calc_avg(void){
        //
        // return the average of whatever score is stored
        //

        vector<Transcript*> vec;

        for (auto itr = this->_transcripts.begin(); itr != this->_transcripts.end(); itr++){
            string transcript_id   = itr->first;
            Transcript &transcript = itr->second;
            transcript.id          = transcript_id; // attach name
            double sum             = 0.0;
            for (uint i = 0; i < transcript.scores.size(); i++){
                sum += transcript.scores[i];
            }
            if (transcript.scores.empty()){
                transcript.avg = 0.0;
            }
            else {
                transcript.avg = sum / (double)transcript.scores.size();
            }
            vec.push_back(&transcript);
        }

        return vec;
    }
    
    bool has_score(void){
        return this->_has_score;
    }

    bool has_coding(void){
        return this->_coding;
    }

    void set_has_coding(void){
        this->_coding = true;
    }
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

void
parse_attributes(string &attributes, vector<Attribute> &atrbVec){

    //
    // parse the attributes string that is expected to
    // be ';' delimited
    //

    if (attributes.empty()){
        return;
    }

    size_t start = 0, next = string::npos, length = attributes.size();
    string part;

    atrbVec.clear(); // ensure new entries

    // iterate over a ';' delimited string
    while (start <= length){
        next = attributes.find(';', start);
        part = next == string::npos ? attributes.substr(start) : attributes.substr(start, next - start);
        
        // strip whitespace
        uint i = 0;
        while (i < part.size() && part[i] == ' '){
            i++;
        }

       part = part.substr(i);

       if (!part.empty()){
            size_t idx = part.find(' '); // find if we have a key value pair
            if (idx != string::npos){
                string key   = part.substr(0, idx);
                string value = part.substr(idx + 1);
                // remove qoutes
                if (value.size() >= 2 && value[0] == '"' && value.back() == '"'){
                    value = value.substr(1, value.size() - 2);
                }
                atrbVec.push_back({key, value});
            }
       }
        // we reached the end
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }
}

int
parse_annotation(const string &ann, unordered_map<string, vector<Gene*>> &genome, 
                 const unordered_map<string, vector<Site>> &markers, 
                 const bool merge_exons, const bool coding){
    //
    // collect only genes on chromosomes with markers
    //

    bool gzipped = is_compressed(ann);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, ann);

    //
    // we will create a mapping
    // for gene_id -> Gene & at
    // the end, move them into 
    // the genome container
    //
    // chr -> gene_id -> Gene
    //
    unordered_map<string, unordered_map<string, Gene*>> gene_map;

    //
    // create the objects we need to store info
    //
    vector<string> parts;
    vector<Attribute> atrbVec;
    Gene *g;
    string line, chrom, gene_id, gene_name, transcript_id;
    uint start, end;
    bool eof = false;

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
        if (line.empty() || line[0] == '#'){
            continue; // skip comment & empty lines
        }
        parse_tabular(line, parts);
        if (parts.size() < 9){
            cerr << "Malformed line in " << ann << '\n' << line;
            exit(1);
        }

        chrom = parts[0];

        if (markers.count(chrom) == 0){
            continue; // no markers on this chrom
        }

        if (parts[2] == "gene" || parts[2] == "exon" || parts[2] == "CDS"){

            if (parts[8].back() == '\n'){
                parts[8].pop_back(); // strip new line char
            }
            // grab the gene_id for this gene
            parse_attributes(parts[8], atrbVec);
            gene_id.clear();
            gene_name.clear();
            transcript_id.clear();
            for (uint i = 0; i < atrbVec.size(); i++){
                if (atrbVec[i].key == "gene_id"){
                    gene_id = atrbVec[i].val;
                }
                else if (atrbVec[i].key == "gene_name"){
                    gene_name = atrbVec[i].val;
                }
                else if (atrbVec[i].key == "transcript_id"){
                    transcript_id = atrbVec[i].val;
                }
            }
            if (gene_id.empty()){
                cerr << "Unable to get gene_id for the following record:\n" << line;
                exit(1);
            }
            chrom = parts[0];
            start = (uint)stoi(parts[3]);
            end   = (uint)stoi(parts[4]);
            if (parts[2] == "gene"){
                g = new Gene(gene_id, gene_name, start, end);
                gene_map[chrom][gene_id] = g;
            }
            else if (gene_map[chrom].find(gene_id) != gene_map[chrom].end()) {
                if (parts[2] == "exon"){
                    Exon exon{start, end};
                    if (transcript_id.empty()){
                        cerr << "Unable to find a transcript_id in the following line:\n"
                             << line << '\n';
                        exit(1);
                    }
                    gene_map[chrom][gene_id]->add_exon(exon, transcript_id);
                }
                else {
                    gene_map[chrom][gene_id]->set_has_coding();
                }
            }
            else {
                cerr << "Malformed annotations. Exon came before gene entry. "
                     << "Offending line:\n" << line;
                exit(1);
            }
        } // end of exon parsing
    } // end of parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (gene_map.empty()){
        cerr << "No genes were found on marker-containing chromosomes in " << ann << '\n';
        exit(1);
    }

    // move the genes into the genome map
    uint total = 0;
    for (auto itr = gene_map.begin(); itr != gene_map.end(); itr++){
        chrom                               = itr->first;
        unordered_map<string, Gene*> &genes = itr->second;
        vector<Gene*> &chrom_genes          = genome[chrom];

        chrom_genes.reserve(genes.size()); // reserve enough space
        for (auto jtr = genes.begin(); jtr != genes.end(); jtr++){

            if (coding && !jtr->second->has_coding()){
                continue; // skip non-coding genes
            }

            if (!merge_exons){
                jtr->second->sort_exons(); // sort exons for each transcript
            }
            else {
                jtr->second->merge_exons(); // collapse gene to a single set of exons
            }
            
            total++;
            chrom_genes.push_back(jtr->second);
        }

        // now sort for downstream binary search
        sort(chrom_genes.begin(), chrom_genes.end(),[]
            (const Gene *geneA, const Gene *geneB){
                if (geneA->start == geneB->start){
                    return geneA->end < geneB->end;
                }
                return geneA->start < geneB->start;
            }  
        );
    }

    cerr << "Loaded " << total << " genes\n";
    return 0;
}

void
get_overlapping_genes(const int left, const int pos, vector<Gene*> *genes, vector<Gene *> &candidates){

    //
    // populate the candidates vector with overlapping genes
    // at this site
    //
    int right = left + 1;

    Gene *gene;

    while (right < genes->size()){
        gene = (*genes)[right];
        if (gene->start <= pos && gene->end >= pos){
            candidates.push_back(gene);
        }
        else if (gene->start > pos){
            break;
        }
        right++;
    }
}

int
parse_pixy(const string &table, const string &pop1, const string &pop2,
            unordered_map<string, vector<Site>> &markers){

    //
    // find sites that are overlapping genes & add each populations
    // score to each gene
    //

    bool gzipped = is_compressed(table);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, table);


    vector<string> parts;
    string line, chrom;
    uint pos1, pos2; // window positions
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

        if (parts.size() < 6){
            cerr << "Error: Expected at least 6 columns. Offending line: " << line << '\n';
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
        pos1  = (uint)stoi(parts[3]);
        pos2  = (uint)stoi(parts[4]);

        // take the middle of the window
        Site site;
        site.pos = (pos2 + pos1) / 2;

        if (parts[5] == "NA" || parts[5][0] == '-'){
            site.value = 0.0;
        }
        else {
            site.value = stod(parts[5]);
        }
      
        markers[chrom].push_back(site);
        found++;

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (found == 0){
        cerr << "Found no markers between " << pop1 << " and " << pop2 << '\n';
        exit(1);
    }

    cerr << "Found " << found << " markers between " << pop1 << " and " << pop2 << '\n';

    return 0;
}

int
overlap_genes(unordered_map<string, vector<Gene*>> &genome,
                   unordered_map<string, vector<Site>> &markers){
    // add each gene to any sites that falls within its boundaries

    uint total = 0;
    for (auto itr = markers.begin(); itr != markers.end(); itr++){
        
        // nothing to overlap here
        if (genome.find(itr->first) == genome.end()){
            continue;
        }
        string chrom         = itr->first;
        vector<Site> &sites  = itr->second;
        vector<Gene*> *genes = &genome[chrom];

        // genes are already sorted, but we need to keep
        // track when a gene ends
        uint furthest = 0, num_genes = genes->size();
        vector<uint> gene_ends(num_genes, 0);
        for (uint i = 0; i < num_genes; i++){
            if (furthest < (*genes)[i]->end){
                furthest = (*genes)[i]->end;
            }
            gene_ends[i] = furthest;
        }

        for (uint i = 0; i < sites.size(); i++){

            // use a binary search to find the leftmost
            // gene this site can overlap
            Gene* gene;
            uint pos  = sites[i].pos;
            uint left = 0, right = num_genes, mid;

            while (left < right){
                mid  = left + (right - left) / 2;
                gene = (*genes)[mid];
                if (gene_ends[mid] >= pos){
                    right = mid;
                }
                else {
                    left = mid + 1;
                }
            }

            bool overlaps = false;
            if (left < num_genes){
                gene     = (*genes)[left];
                overlaps = gene->start <= pos;
            }
            if (overlaps){
                total++;
                vector<Gene *> overlapped_genes = {gene};
                get_overlapping_genes(left, pos, genes, overlapped_genes);
                
                //
                // now add this site's score to each gene
                //
                for (uint j = 0; j < overlapped_genes.size(); j++){
                    overlapped_genes[j]->add_site(sites[i]);
                }
            }

        } // end of sites loop
    } // end of chrom loop

    cerr << "Found " << total << " sites that overlap genes\n";

    return 0;
}

int
write_output(unordered_map<string, vector<Gene *>> &genome, bool const merge_exons){
    //
    // write a tsv where each
    // 
    //

    string gene_id, gene_name;
    string outname = merge_exons ? "Avg_gene_scores.tsv" : "Avg_gene_scores_per_transcripts.tsv";

    ofstream fh;
    fh.open(outname);

    fh << "#GeneID\tGeneName\tTranscriptID\tNumSites\tAvgScore\n";

    // write out all genes, even if they may not carry any heterozygous positions
    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene*> *genes = &itr->second;
        for (uint i = 0; i < genes->size(); i++){
            Gene *gene = (*genes)[i];
            if (!gene->has_score()){
                delete gene;
                continue;
            }
            
            vector<Transcript*> transcripts = gene->calc_avg();
            gene_id   = gene->id;
            gene_name = gene->name.empty() ? gene_id : gene->name;

            fh << setprecision(6);
            for (uint j = 0; j < transcripts.size(); j++){
                fh << gene_id << '\t' << gene_name << '\t'
                   << transcripts[j]->id << '\t'  << transcripts[j]->scores.size() 
                   << '\t' << transcripts[j]->avg << '\n';
            }
            delete gene;        
        }
    }

    fh.close();
    return 0;

}

void
help(){
    cerr << "Usage: ./calc_gene_avg_pixy -t pixy_table -a ann.gtf.gz --pop1 POP1 --pop2 POP2 --merge [optional] --coding [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string table, ann, pop1, pop2;
    bool merge_exons = false, coding = false;
    
    // expect at least 4 inputs
    if (argc < 9){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-t" && i + 1 < argc){
            table = string(argv[i + 1]);
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = string(argv[i + 1]);
        }
        else if (arg == "--pop1" && i + 1 < argc){
            pop1 = string(argv[i + 1]);
        }
        else if (arg == "--pop2" && i + 1 < argc){
            pop2 = string(argv[i + 1]);
        }
        else if (arg == "--merge"){
            merge_exons = true;
        }
        else if (arg == "--coding"){
            coding = true;
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (table.empty() || ann.empty() || pop1.empty() || pop2.empty()){
        help();
    }

    if (!file_exists(table)){
        cerr << "Unable to find " << table << '\n';
        exit(1);
    }
    if (!file_exists(ann)){
        cerr << "Unable to find " << ann << '\n';
        exit(1);
    }

    //
    // first, load all the markers
    //
    unordered_map<string, vector<Site>> markers;
    parse_pixy(table, pop1, pop2, markers);

    //
    // next, load the genes from chromosomes
    // with markers
    //
    unordered_map<string, vector<Gene*>> genome;
    parse_annotation(ann, genome, markers, merge_exons, coding);

    //
    // now overlap the two datasets
    //
    overlap_genes(genome, markers);

    // write the output
    write_output(genome, merge_exons);

    return 0;
}