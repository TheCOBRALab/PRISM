// Iterative HFold files
#include "Result.hh"
#include "W_final.hh"
#include "cmdline.hh"
#include "h_globals.hh"
#include "hotspot.hh"
#include "part_func.hh"
#include "SHAPE.hh"
// a simple driver for the HFold
#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>

bool exists(const std::string path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// check length and if any characters other than ._()
void validateStructure(std::string &seq, std::string &structure) {
    int n = structure.length();
    std::vector<int> pairs;
    for (int j = 0; j < n; ++j) {
        if (structure[j] == '(') pairs.push_back(j);
        if (structure[j] == ')') {
            if (pairs.empty()) {
                std::cerr << "Incorrect input: More left parentheses than right" << std::endl;
                exit(0);
            } else {
                int i = pairs.back();
                pairs.pop_back();
                if (seq[i] == 'A' && seq[j] == 'U') {
                } else if (seq[i] == 'C' && seq[j] == 'G') {
                } else if ((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'G' && seq[j] == 'U')) {
                } else if ((seq[i] == 'U' && seq[j] == 'G') || (seq[i] == 'U' && seq[j] == 'A')) {
                } else if ((seq[i] == 'T' && seq[j] == 'A') || (seq[i] == 'A' && seq[j] == 'T')) {
                } else {
                    std::cerr << "Incorrect input: " << seq[i] << " does not pair with " << seq[j] << std::endl;
                    exit(0);
                }
            }
        }
    }
    if (!pairs.empty()) {
        std::cerr << "Incorrect input: More left parentheses than right" << std::endl;
        exit(0);
    }
}

// check if sequence is valid with regular expression
// check length and if any characters other than GCAUT
void validateSequence(std::string sequence) {

    if (sequence.length() == 0) {
        std::cout << "sequence is missing" << std::endl;
        exit(EXIT_FAILURE);
    }
    // return false if any characters other than GCAUT -- future implement check based on type
    for (char c : sequence) {
        if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T' || c == 'N')) {
            std::cout << "Sequence contains character " << c << " that is not N,G,C,A,U, or T." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

inline void trim(std::string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v\""));
    s.erase(s.find_last_not_of(" \t\n\r\f\v\"") + 1);
}

/**
 * @brief Represents an RNA entry with name, sequence, and structure.
 *
 * This struct is designed to store information about an RNA molecule, 
 * including its name, nucleotide sequence, and secondary structure. 
 * It also provides utility functions for checking the size of the RNA entry.
 */

struct RNAEntry {
    std::string name;
    std::string sequence;
    std::string structure;

    RNAEntry() = default;

    RNAEntry(std::string rna_name, std::string rna_sequence, std::string rna_structure)
        : name{rna_name},
          sequence{rna_sequence},
          structure{rna_structure} {}

    RNAEntry(std::string rna_sequence, std::string rna_structure)
        : RNAEntry("N/A", rna_sequence, rna_structure) {}

    size_t size() const {
        return structure.size();
    }
};

std::vector<RNAEntry> get_all_file_entries(const std::string& file){
    if(!exists(file)){
        std::cerr << "Error: Input file not found: " << file << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // state machine to parse the file
    #define UNINITIALIZED -1
    #define NAME 0
    #define SEQUENCE 1
    #define STRUCTURE 2

    std::ifstream in(file.c_str());
    RNAEntry current;
    std::vector<RNAEntry> entries;
    int state = UNINITIALIZED;
    int line_number = 0;

    std::string str;
    while(getline(in, str)){
        ++line_number;
        trim(str);
        if (str.empty()) continue;

        // Check if the line is the name of the entry
        if ((state == UNINITIALIZED) && (str[0] != '>')) {
            std::cerr << "Error: Expected '>' at the beginning of the line: " << str << ". Line number: " << line_number <<std::endl;
            exit(EXIT_FAILURE);
        }

        if (str[0] == '>'){

            if (!current.name.empty() && !current.sequence.empty() && current.structure.empty()) {
                current.structure = "";
            }

            if (state != UNINITIALIZED) { // valid entry, save it
                if (current.sequence.empty() && current.structure.empty()) {
                    std::cerr << "Warning: Sequence and structure are empty for entry: " << current.name << ". Line number: " << line_number << ". Skipping..."<<  std::endl;
                }
                validateSequence(current.sequence);
                validateStructure(current.sequence,current.structure);
                entries.push_back(current);
                current = RNAEntry();
            }

            current.name = str.substr(1);
            current.sequence = "";
            current.structure = "";
            state = SEQUENCE;

        } else if (state == SEQUENCE){
            if(str.find('(') != std::string::npos || str.find(')') != std::string::npos || str.find('.') != std::string::npos){
                state = STRUCTURE;
            } else{
                current.sequence += str;
            }
            
        } else if (state == STRUCTURE) {
            if(str[0] == '>'){
                state = NAME;
            } else {
                current.structure += str;
            }
        }
    }

    // Handle the last entry
    if (!current.name.empty() && !current.sequence.empty()) {
        if (current.structure.empty()) {
            current.structure = "";
        }
        entries.push_back(current);
    }

    return entries;
}

std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq, const std::string& restricted) {
    std::vector<RNAEntry> entries;
    if (!seq.empty()) {
        entries.emplace_back("Console Sequence", seq, restricted);
    }
    if (!fileI.empty()){
        std::vector<RNAEntry> file_entries = get_all_file_entries(fileI);
        entries.insert(entries.end(), file_entries.begin(), file_entries.end());
    }
    if (entries.empty()) throw std::runtime_error("Sequence is missing");
    return entries;
}

std::string hfold(std::string seq, std::string res, pf_t &energy, sparse_tree &tree, SHAPEData &ShapeData, bool pk_free, bool pk_only, int dangles) {
    W_final min_fold(seq, res,ShapeData, pk_free, pk_only, dangles);
    energy = min_fold.hfold(tree);
    std::string structure = min_fold.structure;
    return structure;
}

std::string hfold_pf(std::string &seq, std::string &final_structure, pf_t &energy, std::string &MEA_structure, pf_t &MEA, std::string &centroid_structure, std::vector<std::pair<std::string,double>> &fatgraphs,pf_t &distance, pf_t &frequency, pf_t &diversity, int &num_fatgraphs, sparse_tree &tree, SHAPEData &ShapeData, bool pk_free,bool pk_only, int dangles, double min_en,
                     int num_samples, bool print_samples, bool PSplot, double gamma) {
    W_final_pf min_fold(seq, final_structure,ShapeData, pk_free,pk_only, dangles, min_en, num_samples,print_samples, PSplot, gamma);
    energy = min_fold.hfold_pf(tree);
    std::string structure = min_fold.structure;
    MEA = min_fold.hfold_MEA(tree);
    MEA_structure = min_fold.MEA_structure;
    distance = min_fold.hfold_centroid(tree);
    centroid_structure = min_fold.centroid_structure;
    min_fold.hfold_fatgraph(fatgraphs,num_fatgraphs);
    diversity = min_fold.ensemble_diversity;
    frequency = min_fold.frequency;
    return structure;
}

void seqtoRNA(std::string &sequence) {
    for (char &c : sequence) {
        if (c == 'T') c = 'U';
    }
}

void print_results(std::vector<Result> &result_list, std::vector<std::vector<std::pair<std::string,double>>> &fatgraphs, std::string &fileO, int number_of_output){
    if (fileO != "") {
        std::ofstream out(fileO,std::fstream::app);
        out << result_list[0].get_sequence() << std::endl;
        for (cand_pos_t i = 0; i < number_of_output; i++) {
            if (i>0 && result_list[i].get_final_structure() == result_list[i - 1].get_final_structure()) continue;
            int fatgraph_num = result_list[i].get_fatgraph_num();
            out << "Restricted_" << i << ": " << result_list[i].get_restricted() << " (" << result_list[i].get_restricted_energy() << ")"
                << std::endl;
            out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")"
                << std::endl;
            out << "Result_" << i << ":     " << result_list[i].get_final_structure_pf() << " (" << result_list[i].get_pf_energy() << ")"
                << std::endl;
            out << "Result_" << i << ":     " << result_list[i].get_MEA_structure() << " (" << result_list[i].get_MEA() << ")" << std::endl;
            out << "Result_" << i << ":     " << result_list[i].get_centroid_structure() << " (" << result_list[i].get_distance() << ")" << std::endl;
            out << "Result_" << i << ":     ";
            for(size_t j=0; j<fatgraphs[fatgraph_num].size();++j){
               out << fatgraphs[fatgraph_num][j].first << "\t(" << fatgraphs[fatgraph_num][j].second << ")\t";
            }
            out << std::endl;
            out << "frequency of MFE structure in ensemble: " << result_list[i].get_frequency() << "; ensemble diversity " << result_list[i].get_diversity() << std::endl;
        }

    } else {
        // kevin: june 22 2017
        // Mateo: Sept 13 2023
        // changed format for ouptut to stdout
        std::cout << result_list[0].get_sequence() << std::endl;
        if (result_list.size() == 1) {
            int fatgraph_num = result_list[0].get_fatgraph_num();
            std::cout << result_list[0].get_restricted() << std::endl;
            std::cout << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
            std::cout << result_list[0].get_final_structure_pf() << " (" << result_list[0].get_pf_energy() << ")" << std::endl;
            std::cout << result_list[0].get_MEA_structure() << " (" << result_list[0].get_MEA() << ")" << std::endl;
            std::cout << result_list[0].get_centroid_structure() << " (" << result_list[0].get_distance() << ")" << std::endl;
            for(size_t j=0; j<fatgraphs[fatgraph_num].size();++j){
               std::cout << std::fixed << std::setprecision(4) << fatgraphs[fatgraph_num][j].first << "\t(" << fatgraphs[fatgraph_num][j].second << ")\t";
            }
            std::cout << std::endl;
            std::cout << "frequency of MFE structure in ensemble: " << result_list[0].get_frequency() << "; ensemble diversity " << result_list[0].get_diversity() << std::endl;
        } else {
            for (cand_pos_t i = 0; i < number_of_output; i++) {
                if (i>0 && result_list[i].get_final_structure() == result_list[i - 1].get_final_structure()) continue;
                int fatgraph_num = result_list[i].get_fatgraph_num();
                std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << " (" << result_list[i].get_restricted_energy() << ")"
                          << std::endl;
                std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")"
                          << std::endl;
                std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure_pf() << " (" << result_list[i].get_pf_energy() << ")"
                          << std::endl;
                std::cout << "Result_" << i << ":     " << result_list[i].get_MEA_structure() << " (" << result_list[i].get_MEA() << ")"
                          << std::endl;
                std::cout << "Result_" << i << ":     " << result_list[i].get_centroid_structure() << " (" << result_list[i].get_distance() << ")"
                          << std::endl;
                std::cout << "Result_" << i << ":     ";
                for(size_t j=0; j<fatgraphs[fatgraph_num].size();++j){
                    std::cout << fatgraphs[fatgraph_num][j].first << "\t(" << fatgraphs[fatgraph_num][j].second << ")\t";
                }
                std::cout << std::endl;
                std::cout << "frequency of MFE structure in ensemble: " << result_list[i].get_frequency() << "; ensemble diversity " << result_list[i].get_diversity() << std::endl;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    args_info args_info;

    // get options (call getopt command line parser)
    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    std::string seq;
    if (args_info.inputs_num > 0) {
        seq = args_info.inputs[0];
    } else {
        if (!args_info.input_file_given) std::getline(std::cin, seq);
    }

    std::string restricted;
    args_info.input_structure_given ? restricted = args_info.input_structure_arg : restricted = "";
    std::string fileI;
    args_info.input_file_given ? fileI = args_info.input_file_arg : fileI = "";

    std::string fileO;
    args_info.output_file_given ? fileO = args_info.output_file_arg : fileO = "";

    int number_of_suboptimal_structure = args_info.subopt_given ? args_info.subopt_arg : 1;

    bool pk_free = args_info.pk_free_given;
    bool pk_only = args_info.pk_only_given;
    std::string shapeFile = args_info.shape_given ? args_info.shape_arg : "";

    int dangles = args_info.dangles_given ? args_info.dangles_arg : 2;

    int num_samples = args_info.samples_given ? args_info.samples_arg : 1000;

    bool print_samples = args_info.print_samples_given;

    double gamma = args_info.gamma_given ? args_info.gamma_arg : 1;

    bool PSplot = !args_info.noPS_given;

    int num_fatgraph = args_info.fatgraph_given ? args_info.fatgraph_arg : 1;

    if(args_info.paramFile_given){
        std::string file = args_info.paramFile_arg;
        if (exists(file)) vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
        else{
            std::cerr << "Not a valid parameter file!" << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        if (seq.find('T') != std::string::npos) {
            vrna_params_load_DNA_Mathews2004();
        } else{
            std::string file = std::string(PARAMS_DIR) + "/rna_DirksPierce09.par";
            if (exists(file)) vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            else{
                std::cerr << "Not a valid parameter file!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    std::vector<RNAEntry> Inputs = get_all_inputs(fileI,seq,restricted);

    for(RNAEntry current : Inputs){
        cand_pos_t n = current.sequence.length();
        std::transform(current.sequence.begin(), current.sequence.end(), current.sequence.begin(), ::toupper);
        if (!args_info.noConv_given) seqtoRNA(current.sequence);
        if (pk_free) if (current.structure == "") current.structure = std::string(n,'.');

        SHAPEData ShapeData(shapeFile,n);

        // Hotspots
        std::vector<Hotspot> hotspot_list;
        vrna_param_s *params;
        params = vrna_params(NULL);
        std::cout << current.sequence << std::endl;
        if (restricted != "") {
            Hotspot hotspot(1, current.structure.length(), current.structure.length() + 1);
            hotspot.set_structure(current.structure);
            hotspot_list.push_back(hotspot);
        }
        if ((number_of_suboptimal_structure - hotspot_list.size()) > 0) {
            get_hotspots(current.sequence, hotspot_list,ShapeData, number_of_suboptimal_structure, params);
        }
        free(params);
        // Data structure for holding the output
        std::vector<Result> result_list;
        //  Iterate through all hotspots or the single given input structure
        cand_pos_t size = hotspot_list.size();
        
        pf_t energy,energy_pf,MEA,distance,frequency,diversity;
        std::string MEA_structure,centroid_structure;
        std::vector<std::vector<std::pair<std::string,double>>> fatgraphs(num_fatgraph);
        for (cand_pos_t i = 0; i < size; ++i) {
            std::string structure = hotspot_list[i].get_structure();
            sparse_tree tree(structure, n);
            std::string final_structure = hfold(current.sequence, structure, energy, tree,ShapeData, pk_free, pk_only, dangles);
            std::string final_structure_pf = hfold_pf(current.sequence, final_structure, energy_pf,MEA_structure,MEA,centroid_structure,fatgraphs[i],distance,frequency, diversity,num_fatgraph, tree,ShapeData, pk_free,pk_only, dangles, energy, num_samples, print_samples, PSplot, gamma);

            if (!args_info.input_structure_given && energy > 0.0) {
                energy = 0.0;
                energy_pf = 0.0;
                final_structure = std::string(n, '.');
            }

            Result result(current.sequence, hotspot_list[i].get_structure(), hotspot_list[i].get_energy(), final_structure, energy, final_structure_pf, energy_pf,MEA_structure,MEA,centroid_structure,distance,i,frequency,diversity);
            result_list.push_back(result);
        }

        Result::Result_comp result_comp;
        std::sort(result_list.begin(), result_list.end(), result_comp);

        int number_of_output = 1;

        if (number_of_suboptimal_structure != 1) {
            number_of_output = std::min((int)result_list.size(), number_of_suboptimal_structure);
        }
        print_results(result_list,fatgraphs,fileO,number_of_output);
    }

    // output to file
    cmdline_parser_free(&args_info);

    return 0;
}
