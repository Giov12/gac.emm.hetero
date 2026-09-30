#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <zlib.h>
#include <iomanip>
#include <algorithm>
#include <unordered_map>

using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
using std::cout;
using std::stoi;
using std::sort;
using std::unordered_map;

//
// this code will take an annotation file
// & create a bed file with the gene boundaries
// & then it will check if a SNP overlaps an exonic
// site. If so, it is written to a sites file
//

typedef unsigned int uint;

struct Exon {
    uint start;
    uint end;
};

struct Attribute {
    string key;
    string val;
};

class Gene {
private:
    bool _coding = false;

public:
    string       id;
    string       name;
    uint         start;
    uint         end;
    vector<Exon> exons;

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
        this->exons.clear();
    }

    void add_exon(Exon exon){
        //
        // just add the exon for now
        // & then we will resolve
        // the exons at the end
        //

        this->exons.push_back(exon);
    }

    void resolve_exons(void){
        //
        // de-duplicate & merge overlapping exons
        //

        // handles empty and single exon cases
        if (this->exons.size() <= 1){
            return;
        }

        vector<Exon> resolved;

        const int count = this->exons.size();
        resolved.reserve(count);

        //
        // sort to then just go exon by exon
        //
        sort(this->exons.begin(), this->exons.end(), []
            (const Exon &exon1, const Exon &exon2){
                if (exon1.start == exon2.start){
                    return exon1.end < exon2.end;
                }
                return exon1.start < exon2.start;
            }
        );

        resolved.push_back(exons.front());
        int i = 1;

        while (i < count){
            Exon &prev = resolved.back();
            Exon &next = this->exons[i];

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
        this->exons = resolved;
    }


    bool is_exonic(uint pos){
        //
        // exons are ordered,
        // so just see if there is an overlap
        //

        for (uint i = 0; i < this->exons.size(); i++){
            const Exon &e = this->exons[i];
            if (e.start <= pos && pos <= e.end){
                return true;
            }
        }

        // no overlap
        return false;
    }

    void set_as_coding(void){
        this->_coding = true;
    }

    bool coding(void) const {
        return this->_coding;
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
make_bed(const string &ann, unordered_map<string, vector<Gene*>> &genome, const bool coding){

    //
    // collect all the genes and also create the bed file in parallel
    //

    cerr << "Parsing " << ann << '\n';

    bool gzipped = is_compressed(ann);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, ann);

    ofstream ofh("genes.bed");

    if (!ofh.is_open()){
        cerr << "ErrorL Unable to write genes.bed in the current directory\n";
        exit(1);
    }

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
    string line, chrom, gene_id, gene_name;
    uint start, end, written = 0;
    bool eof;

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

        if (parts[2] == "gene" || parts[2] == "exon" || parts[2] == "CDS"){

            if (parts[8].back() == '\n'){
                parts[8].pop_back(); // strip new line char
            }
            // grab the gene_id for this gene
            parse_attributes(parts[8], atrbVec);
            gene_id   = "";
            gene_name = "";
            for (uint i = 0; i < atrbVec.size(); i++){
                if (atrbVec[i].key == "gene_id"){
                    gene_id = atrbVec[i].val;
                }
                else if (atrbVec[i].key == "gene_name"){
                    gene_name = atrbVec[i].val;
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
                    gene_map[chrom][gene_id]->add_exon(exon);
                }
                else{
                    gene_map[chrom][gene_id]->set_as_coding();
                }
            }
            else {
                cerr << "Malformed annotations. Exon/CDS came before gene entry. "
                     << "Offending line:\n" << line;
                exit(1);
            }
        } // end of exon parsing
    } // end of parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    if (gene_map.empty()){
        cerr << "No genes were found in " << ann << '\n';
        exit(1);
    }
    // move the genes into the genome map

    for (auto itr = gene_map.begin(); itr != gene_map.end(); itr++){
        chrom                               = itr->first;
        unordered_map<string, Gene*> &genes = itr->second;
        vector<Gene*> &chrom_genes          = genome[chrom];

        chrom_genes.reserve(genes.size()); // reserve enough space
        for (auto jtr = genes.begin(); jtr != genes.end(); jtr++){
            Gene *gene = jtr->second;
            if (coding && !gene->coding()){
                continue; // analysis restricted to protein-coding genes
            }

            // deduplicate & merge overlapping exons
            gene->resolve_exons(); 

            // write the bed file
            string outname = gene->id + '_' +  (gene->name.empty() ? gene->id : gene->name);
            ofh << chrom << '\t' << gene->start - 1 << '\t' << gene->end << '\t' << outname << '\n';
            written++;

            chrom_genes.push_back(gene);
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

    cerr << "Wrote " << written << " gene entries to genes.bed\n";
    return 0;
}

int
make_sites(const string &vcf, unordered_map<string, vector<Gene*>> &genome){

    //
    // this function will be the main work horse for this
    // code
    //

    cerr << "Starting to parse " << vcf << '\n';

    bool gzipped = is_compressed(vcf);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, vcf);

    ofstream ofh("sites.tsv");

    if (!ofh.is_open()){
        cerr << "Error: Unable to open sites.tsv in this directory\n";
        exit(1);
    }
    
    vector<string> parts;
    vector<uint> gene_ends; // for long genes
    vector<Gene*> *genes;
    string chrom, line, curChrom;
    uint pos, ngenes, total = 0, furthest;
    bool  eof = false; // default value
    Gene *gene;

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

        if (line.empty()){
            continue;
        }
        
        if (line[0] == '#'){
            // we do not care about meta data
            continue;
        }


        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }

        parse_tabular(line, parts);

        chrom = parts[0];
        pos   = (uint)stoi(parts[1]);

        // get the new genes & sort out their ending points
        if (chrom != curChrom){
            curChrom = chrom;
            genes    = &genome[chrom];
            ngenes   = genes->size();
            furthest = 0;
            gene_ends.clear();
            gene_ends.resize(ngenes);
            

            for (uint i = 0; i < ngenes; i++){
                if ((*genes)[i]->end > furthest){
                    furthest = (*genes)[i]->end;
                }
                gene_ends[i] = furthest;
            }
        }

        if (ngenes == 0){
            continue; // no genes in this chrom
        }

        // now check if this gene lands in any exonic gene
        uint left = 0, mid, right = ngenes;
        while (left < right){
            mid  = left + (right - left) / 2;
            if (gene_ends[mid] >= pos){
                right = mid;
            }
            else {
                left = mid + 1;
            }
        }
        // did we stop before we got to the end of the list
        bool found = false;
        if (left < ngenes){
            gene  = (*genes)[left]; // gene in the most left boundary
            found = gene->start <= pos;
        }

        if (found){
            // check if it falls within any exons
            found = false; 
            while (left < ngenes){
                gene = (*genes)[left];
                if (gene->start <= pos && pos <= gene->end){
                    if (gene->is_exonic(pos)){
                        found = true;
                        break;
                    }
                }
                else if (pos < gene->start){
                    break;
                }
                left++;
            }
        }

        if (found){
            total++;
            ofh << chrom << '\t' << pos << '\n';
        }

    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);
    ofh.close();

    cerr << "Found a total of " << total << " snps overlapping an exon\n";

    return 0;
}

void
help(){
    cerr << "Usage: ./ann_to_bed_sites -v vcf.gz -a ann.gtf.gz --coding [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, arg;
    bool coding = false;

    // expect at least two arguments
    if (argc < 5){
        help();
    }

    for (int i = 1; i < argc; i++){
        arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = string(argv[i + 1]);
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = string(argv[i + 1]);
        }
        else if (arg == "--coding"){
            coding = true;
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (vcf.empty() || ann.empty()){
        help();
    }

    if (!file_exists(vcf)){
        cerr << "Unable to find " << vcf << '\n';
        exit(1);
    }
    if (!file_exists(ann)){
        cerr << "Unable to find " << ann << '\n';
        exit(1);
    }

    // first, collect the genes
    unordered_map<string, vector<Gene*>> genome;
    make_bed(ann, genome, coding);

    // now find which snps land on exons
    make_sites(vcf, genome);

    return 0;
}