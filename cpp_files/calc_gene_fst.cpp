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

struct SNP {
    uint   pos;
    double fst;
};

struct Attribute {
    string key;
    string val;
};

struct Transcript {
    string id;
    double avg; // average fst
    uint length = 0;
    vector<Exon> exons;
    vector<double> fsts;
};

class Gene {

private:
    // key -> transcript ID
    unordered_map<string, Transcript> _transcripts;
    bool _has_fst = false;

public:
    string id;
    string name;
    uint   start;
    uint   end;

    //
    // empty constructor
    //
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

    void add_fst(SNP &snp){
        //
        // add the fst to each transcript it overlaps
        //
        
        // should not happen
        if (this->_transcripts.empty()){
            return;
        }

        this->_has_fst = true; // overwriting is fine

        for (auto itr = this->_transcripts.begin(); itr != this->_transcripts.end(); itr++){
            vector<Exon> &exons = itr->second.exons;
            for (uint i = 0; i < exons.size(); i++){
                if (exons[i].start <= snp.pos && snp.pos <= exons[i].end){
                    itr->second.fsts.push_back(snp.fst);
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

    vector<Transcript*> calc_fst(void){
        //
        // return the average fst for each transcripts
        //

        vector<Transcript*> vec;

        for (auto itr = this->_transcripts.begin(); itr != this->_transcripts.end(); itr++){
            string transcript_id   = itr->first;
            Transcript &transcript = itr->second;
            transcript.id          = transcript_id; // attach name
            double sum             = 0.0;
            for (uint i = 0; i < transcript.fsts.size(); i++){
                sum += transcript.fsts[i];
            }
            transcript.avg = sum / (double)transcript.length;
            vec.push_back(&transcript);
        }

        return vec;
    }
    
    bool has_fst(void){
        return this->_has_fst;
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
    // get the gene_id from a ';' delimited string
    // if we only want the gene_id
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
    const unordered_map<string, vector<SNP>> &markers, const bool merge_exons){
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

        if (parts[2] == "gene" || parts[2] == "exon"){

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
                Exon exon{start, end};
                if (transcript_id.empty()){
                    cerr << "Unable to find a transcript_id in the following line:\n"
                         << line << '\n';
                    exit(1);
                }
                gene_map[chrom][gene_id]->add_exon(exon, transcript_id);
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

    for (auto itr = gene_map.begin(); itr != gene_map.end(); itr++){
        chrom                               = itr->first;
        unordered_map<string, Gene*> &genes = itr->second;
        vector<Gene*> &chrom_genes          = genome[chrom];

        chrom_genes.reserve(genes.size()); // reserve enough space
        for (auto jtr = genes.begin(); jtr != genes.end(); jtr++){
            if (!merge_exons){
                jtr->second->sort_exons(); // sort exons for each transcript
            }
            else {
                jtr->second->merge_exons(); // collapse gene to a single set of exons
            }
            
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
parse_table(const string &table, unordered_map<string, vector<SNP>> &markers){

    //
    // find snps that are overlapping genes & add each populations
    // score to each gene
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

        if (line_num == 1 && parts.front() == "CHROM"){
            continue; // skip header
        }

        chrom = parts[0];
        SNP snp;
        snp.fst = parts[2][0] == '-' ? 0.0 : stod(parts[2]);
        snp.pos = (uint)stoi(parts[1]);
        markers[chrom].push_back(snp);

        if (snp.fst == 1.0){
            found++;
        }

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    cerr << "Found a total of " << found << " SNPs with an Fst of 1\n";

    return 0;
}

int
overlap_genes(unordered_map<string, vector<Gene*>> &genome,
                   unordered_map<string, vector<SNP>> &markers){
    // add each gene to any SNPs that falls within its boundaries

    uint total = 0;
    for (auto itr = markers.begin(); itr != markers.end(); itr++){
        
        // nothing to overlap here
        if (genome.find(itr->first) == genome.end()){
            continue;
        }
        string chrom         = itr->first;
        vector<SNP> &snps    = itr->second;
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

        for (uint i = 0; i < snps.size(); i++){

            // use a binary search to find the leftmost
            // gene this snp can overlap
            Gene* gene;
            uint pos  = snps[i].pos;
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
                // now add this snp's fst to each gene
                //
                for (uint j = 0; j < overlapped_genes.size(); j++){
                    overlapped_genes[j]->add_fst(snps[i]);
                }
            }

        } // end of snps loop
    } // end of chrom loop

    cerr << "Found " << total << " snps that overlap genes\n";

    return 0;
}

int
write_output(unordered_map<string, vector<Gene *>> &genome, bool const merge_exons){
    //
    // write a tsv where each
    // 
    //

    string gene_id, gene_name;
    string outname = merge_exons ? "Avg_gene_fsts.tsv" : "Avg_gene_fsts_per_transcripts.tsv";

    ofstream fh;
    fh.open(outname);

    fh << "#GeneID\tGeneName\tTranscriptID\tNumSNPS\tAvgFst\n";

    // write out all genes, even if they may not carry any heterozygous positions
    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene*> *genes = &itr->second;
        for (uint i = 0; i < genes->size(); i++){
            Gene *gene = (*genes)[i];
            if (!gene->has_fst()){
                delete gene;
                continue;
            }
            
            vector<Transcript*> transcripts = gene->calc_fst();
            gene_id   = gene->id;
            gene_name = gene->name.empty() ? gene_id : gene->name;

            fh << std::fixed << setprecision(2);
            for (uint j = 0; j < transcripts.size(); j++){
                fh << gene_id << '\t' << gene_name << '\t'
                   << transcripts[j]->id << '\t'  << transcripts[j]->fsts.size() 
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
    cerr << "Usage: ./calc_gene_fst -t fst_table.tsv.gz -a ann.gtf.gz --merge [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string table, ann;
    bool merge_exons = false;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-t" && i + 1 < argc){
            table = string(argv[i + 1]);
        }
        if (arg == "-a" && i + 1 < argc){
            ann = string(argv[i + 1]);
        }
        else if (arg == "--merge"){
            merge_exons = true;
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (table.empty() || ann.empty()){
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
    unordered_map<string, vector<SNP>> markers;
    parse_table(table, markers);

    //
    // next, load the genes from chromosomes
    // with markers
    //
    unordered_map<string, vector<Gene*>> genome;
    parse_annotation(ann, genome, markers, merge_exons);

    //
    // now overlap the two datasets
    //
    overlap_genes(genome, markers);

    // write the output
    write_output(genome, merge_exons);

    return 0;
}