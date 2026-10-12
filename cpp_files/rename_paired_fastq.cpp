#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <cstdint>
#include <vector>
#include <zlib.h>

using std::string;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cerr;
using std::cout;
using std::getline;

//
// code to rename the sequence records for a given
// fastq file from a paired-end sequencing project
// The first file will be used as the namer while
// the second file will be set as the named 
//
//

struct fastq_record {
    string header;
    string seq;
    string spacer;
    string quality;
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
    int idx = path.size() - 1;
    return path[idx - 2] == '.' && path[idx - 1] == 'g' && path[idx] == 'z'; 
}

void 
open_in_filestream(bool gzipped, gzFile &gz_fh, ifstream &fh, const string &infile){
    
    bool bad;
    if (gzipped){
        gz_fh = gzopen(infile.c_str(), "rb");
        bad   = gz_fh == NULL;

        if (!bad){
            gzbuffer(gz_fh, 1 << 20); // 1 Mb buffer
        }
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
close_in_filestream(bool gzipped, gzFile &gz_fh, ifstream &fh, const string &infile){

    bool bad;
    if (gzipped){
        bad = gzclose(gz_fh) != Z_OK;
    }
    else {
        fh.close();
        bad = fh.fail();
    }

    if (bad){
        cerr << "Error: Problem closing " << infile << '\n';
        exit(1);
    }
}

void 
open_out_filestream(bool gzipped, gzFile &gz_fh, ofstream &fh, const string &infile){
    
    bool bad;
    if (gzipped){
        gz_fh = gzopen(infile.c_str(), "wb");
        bad   = gz_fh == NULL;

        if (!bad){
            gzbuffer(gz_fh, 1 << 20);
        }
    }
    else {
        fh.open(infile);
        bad = !fh.is_open();
    }
    if (bad){
        cerr << "Error: could not open " << infile << " for writing\n";
        exit(1);
    }
}

void 
close_out_filestream(bool gzipped, gzFile &gz_fh, ofstream &fh, const string &infile){

    bool bad;
    if (gzipped){
        bad = gzclose(gz_fh) != Z_OK;
    }
    else {
        fh.close();
        bad = fh.fail();
    }

    if (bad){
        cerr << "Error: Unable to finish writing " << infile << ". (disk full?)\n";
        exit(1);
    }
}

string
get_gzline(gzFile &fh, bool &eof){
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

    if (!chars_read){
        // determine if the gzip file was corrupted
        int   error_num = Z_OK;
        const char *msg = gzerror(fh, &error_num);

        if (error_num != Z_OK && error_num != Z_STREAM_END){
            cerr << "Error: gzip read failed " << msg << '\n';
            exit(1);
        }
    }

    eof = !chars_read; // will be true if no characters read
    return line;
}

bool
read_line(bool gzipped, gzFile &gz_fh, ifstream &fh, string &line){
    //
    // return false at the end of the file
    // \n chars are removed
    //

    bool eof;
    if (gzipped){
        line = get_gzline(gz_fh, eof);
        if (eof){
            return false;
        }
    }
    else if (!getline(fh, line)){
        return false;
    }

    // we didn't return, so strip \n
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
        line.pop_back();
    }

    return true;
}

bool
get_record(bool gzipped, gzFile &gz_fh, ifstream &fh, fastq_record &rec, const string &name, uint64_t rec_num){

    //
    // reads in 4 lines (i.e., 1 record) for the current file stream
    // returns false when we reach the end of the file or encounter
    // an error
    //

    if (!read_line(gzipped, gz_fh, fh, rec.header)){
        return false;
    }

    // read the following 3 lines
    if (!read_line(gzipped, gz_fh, fh, rec.seq)    || 
        !read_line(gzipped, gz_fh, fh, rec.spacer) || 
        !read_line(gzipped, gz_fh, fh, rec.quality)){
            cerr << "Error: Record " << name << " ends prematurely (i.e., middle of a record)\n";
            exit(1);
    }

    // ensure the four lines are legal
    if (rec.header.empty() || rec.header[0] != '@'){
        cerr << "Error: Record " << rec_num << " of " << name << " has no '@' leading character in its header\n";
        exit(1);
    }
    if (rec.seq.empty()){
        cerr << "Error: Record " << rec_num << " of " << name << " has an empty sequence\n";
        exit(1);
    }
    if (rec.spacer.empty() || rec.spacer[0] != '+'){
        cerr << "Error: Record " << rec_num << " of " << name << " has no '+' leading character in its spacer\n";
        exit(1);
    }
    if (rec.quality.empty()){
        cerr << "Error: Record " << rec_num << " of " << name << " has no quality scores\n";
        exit(1);
    }
    if (rec.seq.size() != rec.quality.size()){
        cerr << "Error: Record " << rec_num << " of " << name << " has different sequence and quality score lengths\n";
        exit(1);
    }

    return true;
}

int
add_suffix(string &out, const string &header, size_t start, const string &suffix){
    //
    // add the suffix to the header (namer header)
    //
    size_t idx = header.find_first_of(' ', start);

    if (idx == string::npos){
        idx = header.size();
    }

    out.append(header, start, idx - start);
    out += suffix;
    out.append(header, idx, string::npos);

    return 0;
}

int
write_record(bool gzipped, gzFile &gz_fh, ofstream &fh, const string &namer_header,
             const fastq_record &rec, const string &suffix, string &out){
    //
    // write out the record to the output stream
    //

    out.clear(); // remove previous contents & reuse space
    add_suffix(out, namer_header, 0, suffix);
    out += '\n';
    out += rec.seq;
    out += '\n';
    out += rec.spacer;
    out += '\n';
    out += rec.quality;
    out += '\n';

    if (gzipped){
        if (gzwrite(gz_fh, out.data(), (unsigned)out.size()) != (int)out.size()){
            cerr << "Error: write to output failed\n";
            exit(1);
        }
    }
    else {
        fh.write(out.data(), out.size());
        if (!fh){
            cerr << "Error: write to output failed\n";
            exit(1);
        }
    }
    return 0;
}
int
rename_fastq(const string &fastq1, const string &fastq2, const string &outfile, const string &suffix){

    //
    // function to use fastq1 sequence reads as the naming records
    // to reheader the records in fastq2 and place those records into the
    // output
    //

    bool f1_gz = is_compressed(fastq1);
    bool f2_gz = is_compressed(fastq2);
    bool of_gz = is_compressed(outfile);
    gzFile f1_gz_fh = NULL, f2_gz_fh = NULL, out_gz_fh = NULL;
    ifstream f1_tx_fh, f2_tx_fh;
    ofstream out_tx_fh;

    open_in_filestream(f1_gz, f1_gz_fh, f1_tx_fh, fastq1);
    open_in_filestream(f2_gz, f2_gz_fh, f2_tx_fh, fastq2);
    open_out_filestream(of_gz, out_gz_fh, out_tx_fh, outfile);


    //
    // instead of reading all the headers for fastq1 into memory
    // we will incrementally move down the files, renaming as we go
    //
    string out; // out is a buffer
    fastq_record rec_f1, rec_f2;
    uint64_t n = 0;
    bool keep_reading_f1, keep_reading_f2;

    // reserve some space
    out.reserve(1024);

    while (true){
        keep_reading_f1 = get_record(f1_gz, f1_gz_fh, f1_tx_fh, rec_f1, fastq1, n + 1);
        keep_reading_f2 = get_record(f2_gz, f2_gz_fh, f2_tx_fh, rec_f2, fastq2, n + 1);

        if (!keep_reading_f1 && !keep_reading_f2){
            break; // we are done reading both files
        }
        if (keep_reading_f1 != keep_reading_f2){
            string early, late;
            if (keep_reading_f1){
                early = fastq2;
                late  = fastq1;
            }
            else {
                early = fastq1;
                late  = fastq2;
            }
            cerr << "Error: " << early << " ended before " << late << '\n';
            exit(1);
        }

        n++;
        write_record(of_gz, out_gz_fh, out_tx_fh, rec_f1.header, rec_f2, suffix, out);
    }

    close_in_filestream(f1_gz, f1_gz_fh, f1_tx_fh, fastq1);
    close_in_filestream(f2_gz, f2_gz_fh, f2_tx_fh, fastq2);
    close_out_filestream(of_gz, out_gz_fh, out_tx_fh, outfile);
    
    cerr << "Renamed " << n << " records into " << outfile << '\n';

    return 0;
}

void
help(){
    cerr << "Usage: ./rename_paired_fastq -1 namer.fastq[.gz] -2 named.fastq[.gz] -o out.fastq[.gz] -s suffix [default = /2]\n";
    cerr << "Fastq files must have the same number of records\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string fastq1, fastq2, outfile, arg;
    string suffix = "/2";
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        arg = argv[i];
        if (arg == "-1" && i + 1 < argc){
            fastq1 = argv[i + 1];
        }
        else if (arg == "-2" && i + 1 < argc){
            fastq2 = argv[i + 1];
        }
        else if (arg == "-o" && i + 1 < argc){
            outfile = argv[i + 1];
        }
        else if (arg == "-s" && i + 1 < argc){
            suffix = argv[i + 1];
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (fastq1.empty() && fastq2.empty() && outfile.empty() && suffix.empty()){
        help();
    }
    if (fastq1.empty()){
        cerr << "-1 is required\n";
        exit(1);
    }
    else if (!file_exists(fastq1)){
        cerr << "Could not find " << fastq1 << '\n';
        exit(1);
    }
    if (fastq2.empty()){
        cerr << "-2 is required\n";
        exit(1);
    }
    else if (!file_exists(fastq2)){
        cerr << "Could not find " << fastq2 << '\n';
        exit(1);
    }
    if (outfile.empty()){
        cerr << "Please specify name of output file\n";
        exit(1);
    }

    //
    // a single function will orchestrate everything
    //
    rename_fastq(fastq1, fastq2, outfile, suffix);
    
    return 0;
}