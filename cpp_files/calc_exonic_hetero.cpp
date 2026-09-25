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
using std::sort;
using std::find;
using std::setprecision;

typedef unsigned int uint;

struct Exon {
    uint start;
    uint end;
};

class Gene {
public:
    string       id;
    uint         start;
    uint         end;
    vector<Exon> exons;
    unordered_map<string, vector<double>> sites;

    //
    // empty constructor
    //
    Gene (string id_, uint start, uint end){
        this->id    = id_;
        this->start = start;
        this->end   = end;
    };
    ~Gene(){
        this->exons.clear();
        this->sites.clear();
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

        // edge-case
        if (this->exons.empty()){
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

    void add_pop_genos(const uint pos, const string pop, const double percent_hetero){

        //
        // tally the number of heteros for this population
        // if the position falls within an exon
        //
        bool valid = false;

        for (uint i = 0; i < this->exons.size(); i++){
            if (this->exons[i].start <= pos && this->exons[i].end >= pos){
                valid = true;
                break;
            }
        }

        // must be intronic
        if (!valid){
            return;
        }

        this->sites[pop].push_back(percent_hetero);
    }

    unordered_map<string, double> calc_hetero(void) const {
        //
        // return the level of heterozygosity at this gene
        // for each pop
        //
        unordered_map<string, double> gene_hetero;

        for (auto itr = this->sites.begin(); itr != this->sites.end(); itr++){
            const string &pop              = itr->first;
            const vector<double> &hets_vec = itr->second;

            if (hets_vec.empty()){
                gene_hetero[pop] = -1.0; // nothing called
                continue;
            }
            double het_sum = 0;
            double total   = hets_vec.size();

            for (uint i = 0; i < hets_vec.size(); i++){
                het_sum += hets_vec[i];
            }
            gene_hetero[pop] = het_sum / total;
        }

        return gene_hetero;
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

string
get_geneid(string &attributes){

    //
    // get the gene_id from a ';' delimited string
    //

    if (attributes.empty()){
        return "";
    }

    size_t start = 0, next = string::npos, length = attributes.size();
    string part;

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
                if (key == "gene_id"){
                    return value; // found it
                }
            }
       }
        // we reached the end
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }

    return "";
}

int
parse_annotation(const string &ann, unordered_map<string, vector<Gene*>> &genome){
    //
    // collect all the genes and population the genome
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
    Gene *g;
    string line, chrom, gene_id;
    uint start, end;
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
        if (parts[2] == "gene" || parts[2] == "exon"){
            // grab the gene_id for this gene
            gene_id = get_geneid(parts[8]);
            if (gene_id.empty()){
                cerr << "Unable to get gene_id for the following record:\n" << line;
                exit(1);
            }
            chrom = parts[0];
            start = (uint)stoi(parts[3]);
            end   = (uint)stoi(parts[4]);
            if (parts[2] == "gene"){
                g = new Gene(gene_id, start, end);
                gene_map[chrom][gene_id] = g;
            }
            else {
                Exon exon{start, end};
                gene_map[chrom][gene_id]->add_exon(exon);
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
            jtr->second->resolve_exons(); // deduplicate & collapse overlapping exons
            chrom_genes.push_back(std::move(jtr->second));
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

double
calc_percent_hetero(string &counts){

    //
    // parse a '/' delimited string
    //

    int homo_ref = 0;
    int hetero   = 0;
    int homo_alt = 0;
    int start    = 0;
    int end      = 0;
    int cnt      = 0;
    int tot      = 0;
    int val;

    while (end < counts.size()){
        if (counts[end] == '/'){
            val = stoi(counts.substr(start, end - start));
            cnt++;
            
            switch (cnt)
                {
                case 1:
                    homo_ref = val;
                    break;
                case 2:
                    hetero = val;
                    break;
                case 3:
                    homo_alt = val; // should not happen
                    break;
                default:
                    break;
                }
            start = end + 1;
        }
        end++;
    }

    if (start < counts.size()){
        val = stoi(counts.substr(start, end - start));
        if (cnt != 2){
            cerr << "Invalid entry found: " << counts << '\n';
            exit(1);
        }
        homo_alt = val;
    }

    tot = homo_ref + homo_alt + hetero;

    // no genotype information for this population
    if (tot == 0){
        return -1.0;
    }

    return hetero / tot;

}

void
get_overlapping_genes(const int middle, const int site, vector<Gene*> *genes, vector<Gene *> &candidates){

    //
    // populate the candidates vector with overlapping genes
    // at this site
    //
    int right = middle + 1;

    Gene *gene;

    while (right < genes->size()){
        gene = (*genes)[right];
        if (gene->start <= site && gene->end >= site){
            candidates.push_back(gene);
        }
        else if (gene->start > site){
            break;
        }
        right++;
    }
}

int
parse_table(const string &table, unordered_map<string, vector<Gene*>> &genome,
            vector<string> &pops){

    //
    // find snps that are overlapping genes & add each populations
    // score to each gene
    //

    bool gzipped = is_compressed(table);
    gzFile gz_fh = NULL;
    ifstream txt_fh;

    open_in_filestream(gzipped, gz_fh, txt_fh, table);


    vector<string> parts;
    vector<Gene*> *genes;
    vector<uint> gene_ends; // to handle long genes
    Gene *gene;
    string line, counts, pop, chrom, curChrom;
    bool   eof;

    //
    // counters to find the column &
    // number of variant sites within exons
    //
    long overlapping = 0, line_num = 0;

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

        parse_tabular(line, parts);

        if (line_num == 1 && parts.size() > 2 && parts.front() == "Chr"){
            for (uint i = 2; i < parts.size(); i++){
                pops.emplace_back(parts[i]);
            }
            continue;
        }
        else {
            cerr << "Expected header: Chr\tPos\tPop1\tPop2.. in " << table << '\n';
            exit(1);
        }

        if (line.empty()){
            continue;
        }
        parse_tabular(line, parts);

        chrom = parts[0];

        if (curChrom != chrom){
            // will trigger at the beginning
            // or at every transition
            curChrom = chrom;
            genes    = &genome[chrom];

            // collect the gene ends
            gene_ends.clear();
            gene_ends.resize(genes->size());
            uint longest = 0;
            for (uint i = 0; i < genes->size(); i++){
                longest = longest > (*genes)[i]->end ? 
                          longest : (*genes)[i]->end;
                gene_ends[i] = longest;
            }
        }
        if (genes->empty()){
            continue; // no genes
        }

        // binary search for all genes that this sites can be located on
        uint site = stoi(parts[1]);
        int left  = 0, mid;
        int right = genes->size();
        while (left < right){
            mid  = left + (right - left) / 2;
            gene = (*genes)[mid];
            if (gene_ends[mid] >= site){
                right = mid;
            }
            else {
                left = mid + 1;
            }
        }

        bool found = false;
        if (left < genes->size()){
            // we stopped somewhere
            gene  = (*genes)[left];
            found = gene->start <= site;
            mid   = left;
        }
        if (found){
            overlapping++;

            // we need to get the range of genes this site covers
            vector<Gene *> within_range;

            within_range.push_back((*genes)[mid]);
            get_overlapping_genes(mid, site, genes, within_range);

            // now we need to supply each gene
            // with the genotypes for each population
            double percnt;
            for (uint i = 2; i < parts.size(); i++){
                pop    = pops[i - 2];
                counts = parts[i];
                percnt = calc_percent_hetero(counts);

                if (percnt == -1){
                    continue; // no genotype info
                }
                for (uint j = 0; j < within_range.size(); j++){
                    within_range[j]->add_pop_genos(site, pop, percnt);
                }
            }   

        }
    } // end of file parsing

    close_in_filestream(gzipped, gz_fh, txt_fh);

    cerr << "Found a total of " << overlapping << " sites overlapping one or more genes\n";

    return 0;
}

int
write_output(unordered_map<string, vector<Gene *>> &genome, vector<string> &pops){
    //
    // write a tsv where each
    // population will have its estimated
    // heterozygosity at each gene written
    //

    unordered_map<string, double> gene_hetero;

    ofstream fh;
    fh.open("Gene.Heterozygosity.tsv");

    fh << "#Gene";
    for (uint j = 0; j < pops.size(); j++){
        fh << '\t' << pops[j];
    }
    fh << '\n';

    // write out all genes, even if they may not carry any heterozygous positions
    for (auto itr = genome.begin(); itr != genome.end(); itr++){
        vector<Gene*> *genes = &itr->second;
        for (uint i = 0; i < genes->size(); i++){
            Gene *gene = (*genes)[i];
            fh << gene->id;

            //
            // get each populations heterozygosity for this gene
            //
            gene_hetero = gene->calc_hetero();
            
            // now write it out to the table
            cerr << std::fixed << std::setprecision(5);
            for (uint k = 0; k < pops.size(); k++){
                auto ktr      = gene_hetero.find(pops[k]); // when no sites found
                double hetero = ktr != gene_hetero.end() ? ktr->second : -1.0;
                fh << '\t' << hetero; 
            }
            fh << '\n';
            // we are done with this gene
            delete gene;            
        }
    }

    fh.close();
    return 0;

}

void
help(){
    cerr << "Usage: ./calc_exonic_hetero -t merged_hwe.tsv.gz -a ann.gtf.gz\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string table, ann;
    
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
    // load the genes
    //
    unordered_map<string, vector<Gene*>> genome;
    parse_annotation(ann, genome);

    // parse the table
    vector<string> pops;
    parse_table(table, genome, pops);

    // write the output
    write_output(genome, pops);

    return 0;
}