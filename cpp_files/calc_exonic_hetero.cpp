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

    void add_pop_genos(const uint pos, const string pop, vector<string> &genos){

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

        if (genos.empty()){
            // due to missing data
            this->sites[pop].push_back(-1.0);
            return;
        }

        uint hets = 0; 
        for (uint i = 0; i < genos.size(); i++){
            const string &geno = genos[i];
            // either delimiter is valid
            size_t p = geno.find_first_of("/|");
            if (p == string::npos){
                continue; // not valid, haploid call?
            }
            if (geno.substr(0, p) != geno.substr(p + 1)){
                hets++; // different alleles
            }
        } // end of i

        // get fraction of hetero individuals
        this->sites[pop].push_back(((double)hets / genos.size()));
    }

    unordered_map<string, double> calc_hetero(void) const{
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
            double total = 0;
            uint  called = 0;

            for (uint i = 0; i < hets_vec.size(); i++){
                if (hets_vec[i] == -1.0){
                    continue; // missing
                }
                total += hets_vec[i];
                called++;
            }
            gene_hetero[pop] = called == 0 ? -1 : total / (double)called;
        }

        return gene_hetero;
    }  
};

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
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

    gzFile fh = gzopen(ann.c_str(), "rb");

    if (fh == NULL){
        cerr << "Error: could not open " << ann << '\n';
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
    Gene *g;
    string line, chrom, gene_id;
    uint start, end;
    bool eof;

    while (true){
        line = get_gzline(fh, eof);

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

    gzclose(fh);

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

int
create_popmap(const string &popmap_file, unordered_map<string, string> &popmap){
    //
    // read a simple tsv of the population map to assign individuals to
    //
    fstream fh(popmap_file);

    if (!fh.is_open()){
        cerr << "Error: Unable to open " << popmap_file << '\n';
        exit(1);
    }

    vector<string> parts, seen;
    string pop, sample, line;
    while (std::getline(fh, line)){
        if (line.empty() || line[0] == '#'){
            continue;
        }
        if (line.back() == '\n'){
            line.pop_back();
        }
        parse_tabular(line, parts);
        if (parts.size() != 2){
            cerr << "Error: Unexpected line in population map: " << popmap_file << '\n';
            cerr << line << '\n';
            exit(1);
        }
        pop            = parts[0];
        sample         = parts[1];
        popmap[sample] = pop;
        if (find(seen.begin(), seen.end(), pop) == seen.end()){
            seen.push_back(pop);
        }
    }

    fh.close();

    cerr << "Loaded " << popmap.size() << " samples across " << seen.size() << " populations\n";

    return 0;
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
parse_vcf(const string &vcf, unordered_map<string, string> &popmap,
          unordered_map<string, vector<Gene*>> &genome){

    gzFile fh = gzopen(vcf.c_str(), "rb");
    if (fh == NULL){
        cerr << "Error: could not open " << vcf << '\n';
        exit(1);
    }

    //
    // population -> sample_idx, sample_idx
    // where sample_idx is the column index
    // for that sample's genotype in the vcf
    //
    unordered_map<string, vector<uint>> pop_indices;

    vector<string> parts;
    vector<Gene*> *genes;
    vector<uint> gene_ends; // to handle long genes
    Gene *gene;
    string line, sample, pop, chrom, curChrom;
    bool   eof;

    //
    // calculate heterozygosity at every exonic
    // site for each population
    //
    long overlapping = 0;

    while (true){
        line = get_gzline(fh, eof);

        if (eof){
            break;
        }

        if (line[0] == '#'){
            //
            // we may need to collect the sample indices
            // check for column starting with #CHROM 
            //
            if (line.size() > 6 && line.substr(0, 6) == "#CHROM"){
                if (line.back() == '\n'){
                    line.pop_back(); // we need to parse this file
                    parse_tabular(line, parts);
                    if (parts.size() < 10){
                        // there is no genotype info here
                        cerr << "No sample information found in " << vcf << '\n';
                        exit(1);
                    }
                    // now parse through
                    for (uint i = 9; i < parts.size(); i++){
                        sample = parts[i];
                        //
                        // this may cause problems
                        // in the future, it is not flexible
                        //
                        pop = popmap[sample];
                        pop_indices[pop].push_back(i);
                    }
                }
            }
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
            for (auto itr = pop_indices.begin(); itr != pop_indices.end(); itr++){
                pop                   = itr->first;
                vector<uint> &indices = itr->second;
                vector<string> genos;
                for (uint i = 0; i < indices.size(); i++){
                    string geno = parts[indices[i]];
                    geno        = geno.substr(0, geno.find(':')); 
                    if (geno.find('.') == string::npos){ // skip missing data?
                        genos.push_back(geno);
                    }
                }
                for (uint i = 0 ; i < within_range.size(); i++){
                    within_range[i]->add_pop_genos(site, pop, genos);
                }
            }

        }
    } // end of file parsing

    gzclose(fh);

    cerr << "Found a total of " << overlapping << " sites overlapping one or more genes\n";

    return 0;
}

int
write_output(unordered_map<string, vector<Gene *>> &genome, unordered_map<string, string> &popmap){
    //
    // write a tsv where each
    // population will have its estimated
    // heterozygosity at each gene written
    //

    unordered_map<string, double> gene_hetero;

    // create a header
    vector<string> pops;
    for (auto jtr = popmap.begin(); jtr != popmap.end(); jtr++){
        if (find(pops.begin(), pops.end(), jtr->second) == pops.end()){
            pops.push_back(jtr->second);
        }
    }

    //
    // ensure consistency
    //
    sort(pops.begin(), pops.end());
    
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
    cerr << "Usage: ./calc_exonic_hetero -v vcf.gz -a annotation.gtf.gz -p popmap.tsv\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string vcf, ann, popmap_file;
    
    // expect at least 3 inputs
    if (argc < 4){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-v" && i + 1 < argc){
            vcf = string(argv[i + 1]);
        }
        else if (arg == "-a" && i + 1 < argc){
            ann = string(argv[i + 1]);
        }
        else if (arg == "-p" && i + 1 < argc){
            popmap_file = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (vcf.empty() && ann.empty() && popmap_file.empty()){
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
    if (!file_exists(popmap_file)){
        cerr << "Unable to find " << popmap_file << '\n';
        exit(1);
    }

    //
    // create the popmap
    //
    unordered_map<string, string> popmap;
    create_popmap(popmap_file, popmap);

    //
    // load the genes
    //
    unordered_map<string, vector<Gene*>> genome;
    parse_annotation(ann, genome);

    // parse the vcf file
    parse_vcf(vcf, popmap, genome);

    // write the output
    write_output(genome, popmap);

    return 0;
}