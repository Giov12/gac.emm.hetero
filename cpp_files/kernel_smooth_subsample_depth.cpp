#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <math.h>
#include <algorithm>

//
// the methodology of this approach was adopted from Catchen et al. (2013)
// Catchen, J., Hohenlohe, P. A., Bassham, S., Amores, A., & Cresko, W. A. (2013).
// Stacks: an analysis tool set for population genomics. Molecular ecology, 22(11), 3124-3140.
//

/*
Had to compile using the following flags -lz
*/

using std::string;
using std::vector;
using std::cout;
using std::stoi;

// sigma & weights will be global
const int SIGMA = 150000;
const int window_span = (3 * SIGMA) + 1;
double weights[window_span];

// create a global file stream
std::ofstream ofh;

typedef unsigned int uint;

struct Record {
    string chrom;
    uint   pos;
    uint   depth;
};

void 
calc_weights(double (&weights)[window_span]){

    // fill weights once
    for (int i = 0; i < window_span; i++){
        weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(SIGMA, 2)));
    }
}

int 
median(vector<uint> &temp){
    //
    // return the median value from this vector
    //
    int median, mid;
    int arraysize = temp.size();

    if (arraysize == 1){
        return temp[0]; }

    mid = arraysize / 2;

    //sort vector
    std::sort(temp.begin(), temp.end());

    if (arraysize % 2 == 0){ // even
        median = (temp[mid] + temp[mid - 1]) / 2;
    }
    else{ // odd
        median = temp[mid];
    }

    return median;
}

void
parse_tabular(const string &line, Record &rec){
    //
    // file a record with the values from this line
    //

    int start = 0, end = 0, length = line.size();

     while (end < length){
        if (line[end] == '\t'){
            if (rec.chrom.empty()){
                rec.chrom = line.substr(start, end - start);
            }
            else {
                rec.pos = stoi(line.substr(start, end - start));
            }
            start = end + 1;
        }
        else if (line[end] == '\n'){
            rec.depth = stoi(line.substr(start, end - start));
        }
        end++;
    }

}

void 
create_smooth_array(vector<uint> &positions, vector<uint> &depths, string chrom){
    
    /*
    positions      = bp positions
    depths         = depths at each position
    weights        = precomputed kernel smoothed weights
    start_pos      = beginning of window on chrom
    end_pos        = end of window on chrom
    window_cen     = current position on chrom
    window_start   = used to find start_pos by adjusting from center
    window_end     = used to find end_pos by adjusting to window_start
    n_window_sites = number of positions
    distnace       = distance from window center
    pos_val        = depth at position
    pos_num        = position in bp
    window_avg     = will hold weighted average
    window_weight  = will accumulate weights to use for weighing average
    weight         = kernel weight at position
    */ 

    int start_pos = 0, end_pos = 0, n_sites = positions.size();
    int window_cen, window_start, window_end, n_window_sites, distance, pos_val, pos_num;
    double window_avg, window_weight, weight;

    vector<double> smth_ar(n_sites);

    for (int site = 0; site < n_sites; site++){
        window_cen   = positions[site]; // go to center
        window_start = 1 <= window_cen - SIGMA ? window_cen - SIGMA : 1; // go to start
        if (window_start < 0){ 
            window_start = 1; //keep within bound
        }
        window_end = window_cen + SIGMA; // go to end
        while (positions[start_pos] < window_start){ // start pointer in positions
            start_pos++;
        }
        while (end_pos < n_sites && positions[end_pos] < window_end){ // end pointer in positions
            end_pos++;
        }

        n_window_sites = 0;   // window size
        window_avg     = 0.0; // reset averages
        window_weight  = 0.0;
        
        for (int pos = start_pos; pos < end_pos; pos++){
            n_window_sites++;
            pos_val  = depths[pos];                    // get depth
            pos_num  = positions[pos];                 // get current position
            distance = (int)abs(pos_num - window_cen); // relative to center
            
            if (distance >= window_span){
                continue; }
            
            weight         = weights[distance];  // get weight at position
            window_avg    += (pos_val * weight); // accumulate sum
            window_weight += weight;             // accumulate weight
        }
        
        if (n_window_sites == 0){
            exit(EXIT_FAILURE);
        }
        
        window_avg    = window_avg/window_weight; // adjust
        smth_ar[site] = window_avg; // set

    }

    // now to write to output
    for (int i = 0; i < positions.size(); i++){
        ofh << chrom << '\t' << positions[i] << '\t' << smth_ar[i] << '\n'; }

}

bool 
is_compressed(const string file_name){
    //
    // all this function is doing is check if the file
    // ends with gz
    //

    if (file_name.size() < 4){
        return false;
    }

    int pos = file_name.size() - 1;

    if (file_name[pos] != 'z' || file_name[pos - 1] != 'g' || file_name[pos - 2] != '.'){
        return false;
    }

    return true;
}

void 
create_filestream(const char *depth_file){

    string outname;
    string fileName = depth_file;

    // get basename of file
    string file_basename =
        fileName.substr(fileName.find_last_of("/\\") + 1);

    // add suffix
    outname = "Smoothed_" + file_basename.substr(0, file_basename.find_last_of('.'));

    ofh = std::ofstream(outname, std::ios::out);

}

void
parse_compressed(const char *depth_file, int subsample){

    int buffer_line_size = 80;                   // read in size
    struct gzFile_s *gzf;                        // pointer to file
    char line_buffer[buffer_line_size];          // char array to hold line
    gzf      = gzopen(depth_file, "r");          // file handle
    bool eol = false;                            // end of file
    string line;                                 // char arr -> string
    memset(line_buffer, '\0', buffer_line_size); // fill in nulls
    string cur_chrom;
    uint median_depth;
    vector<uint> positions;
    vector<uint> depths;
    vector<uint> temp_depth; // hold temp depths
    vector<uint> temp_pos; // hold temp pos

    // fill weights
    calc_weights(weights);

    do
    {
        if (gzgets(gzf, line_buffer, buffer_line_size) == NULL) {
            break;
        }
        if (line_buffer[0] == '\n'){
            continue;
        }
        else {
            // convert char array to string
            line = string(line_buffer);
            Record rec;

            // parse the line
            parse_tabular(line, rec);

            if (cur_chrom.empty()){
                cur_chrom = rec.chrom; }

            if (cur_chrom == rec.chrom){
                temp_depth.push_back(rec.depth);
                temp_pos.push_back(rec.pos);

                if (temp_depth.size() == subsample){
                    median_depth = median(temp_depth);
                    depths.push_back(median_depth);
                    positions.push_back(temp_pos.back());
                    temp_depth.clear();
                    temp_pos.clear();
                }
            }
            
            else{
                if (!temp_depth.empty()){
                    median_depth = median(temp_depth);
                    depths.push_back(median_depth);
                    positions.push_back(temp_pos.back());
                }
                if (!positions.empty()){
                    create_smooth_array(positions, depths, cur_chrom); }
                depths.clear(); // clear objects
                positions.clear();
                temp_depth.clear();
                temp_pos.clear();
                cur_chrom = rec.chrom;
                temp_depth.push_back(rec.depth);
                temp_pos.push_back(rec.pos);
            }     
            } 
    } while (!gzeof(gzf) && !eol);

    // get last chrom
    if (temp_depth.size() > 0){
        median_depth = median(temp_depth);
        depths.push_back(median_depth);
        positions.push_back(temp_pos.back()); // last called position
        create_smooth_array(positions, depths, cur_chrom);
        depths.clear();
        positions.clear();
        temp_depth.clear();
        temp_pos.clear();
    }
}

int 
main(int argc, char *args[]){
    // params
    string depth_file;
    string output_file;
    string message;
    string file_basename;
    int subsample = 1000; // default

    // in case of help or error
    message = "./kernel_smooth -d depth_file [required] -s int [optional]\n";

    // read args
    string param;

    for (int i = 0; i < argc; i++) {
        param = string(args[i]);

        if (param == "-h" || param == "--help") {
            cout << message;
            return 0;
        }
        else if (param == "-d" || param == "--depth_file"){
            depth_file = string(args[i + 1]);
        }
        else if (param == "-s" || param == "--subsample"){
            subsample = std::stoi(args[i + 1]);
        }
    }

    if (depth_file.empty()){
        cout << "No depth file supplied.\n";
        cout << message;
        exit(1);
    }

    // check if gzipped b/c that's how i wrote the program
    if (is_compressed(depth_file)){
        // process
        const char *dfile = depth_file.c_str();
        create_filestream(dfile);
        parse_compressed(dfile, subsample);
    }
    else{
        cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
        exit(1);
    }
    
    return 0;
}