#include "base_types.hh"
#include "part_func.hh"

#include <string>
#include <iostream>
#include <vector>


/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
std::string W_final_pf::compute_centroid(sparse_tree &tree, pf_t &dist, pf_t &diversity){
    dist = 0;
    diversity = 0;
    pf_t p = 0;
    std::string centroid = std::string(n, '.');

    for (cand_pos_t i = 1; i <= n; i++){
        for (cand_pos_t j = i + 1; j <= n; j++) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            p = (pf_t)samples[base_pair] / num_samples;
            diversity += p*(1.0-p);
            if (p > 0.5) {
                /* regular base pair */
                if(tree.weakly_closed(i,j)){
                    centroid[i - 1] = '(';
                    centroid[j - 1] = ')';
                }
                else{
                    centroid[i - 1] = '[';
                    centroid[j - 1] = ']';
                }
                dist += (1 - p);
            } else {
                dist += p;
            }
        }
    }
    diversity*=2; // As there are two sides of a base pair
    return centroid;
}


/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
std::string W_final_pf::compute_centroid_PK_only(sparse_tree &tree, pf_t &dist, pf_t &diversity){
    dist = 0;
    diversity = 0;
    pf_t p = 0;
    std::string centroid = std::string(n, '.');

    // I have to recalculate sample frequencies using only structures that contain pseudoknots
    std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples_PK;
    int num_samples_PK = 0;
    for(auto &it: structures){
        std::string structure = it.first;
        int frequency = it.second;
        if((structure.find('[') != std::string::npos)){
            cand_pos_t length = structure.length();
            num_samples_PK += frequency;
            std::vector<int> paren;
            std::vector<int> sb;
            for(cand_pos_t j=0;j<length;++j){
                if(structure[j] == '(') {
                    paren.push_back(j);
                    continue;
                }
                if(structure[j] == '[') {
                    sb.push_back(j);
                    continue;
                }

                if (structure[j] == ')'){
                    int x = paren[paren.size()-1];
                    paren.pop_back();
                    std::pair<cand_pos_tu, cand_pos_tu> base_pair(x+1,j+1);
                    samples_PK[base_pair]+=frequency;
                }
                if (structure[j] == ']'){
                    int x = sb[sb.size()-1];
                    sb.pop_back();
                    std::pair<cand_pos_tu, cand_pos_tu> base_pair(x+1,j+1);
                    samples_PK[base_pair]+=frequency;
                }
            }
        }
    }

    //Calculate centroid based on PK samples
    for (cand_pos_t i = 1; i <= n; i++){
        for (cand_pos_t j = i + 1; j <= n; j++) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            p = (pf_t)samples_PK[base_pair] / num_samples_PK;
            diversity += p*(1.0-p);
            if (p > 0.5) {
                /* regular base pair */
                if(tree.weakly_closed(i,j)){
                    centroid[i - 1] = '(';
                    centroid[j - 1] = ')';
                }
                else{
                    centroid[i - 1] = '[';
                    centroid[j - 1] = ']';
                }
                dist += (1 - p);
            } else {
                dist += p;
            }
        }
    }
    diversity*=2; // As there are two sides of a base pair
    return centroid;
}
/**
 * If I have an unpaired matrix, I can try to walk through one and get all pairings.
 * The second pass, I will do the same stack for pushing, but I will also have another for the fatgraph
 * If either stack is empty and you get to and opening base (op), you push it into the string. If the distance between the two bases of the same
 * type is filled with only unpaired relative to its parent, you skip adding it to the fatgraph. Popping from the stack does not affect this calculation
 * because as long as the stack has something inside that an internal could form from, popping won't affect it.
 * So when adding another op to the stack when the stack is not empty. I get the cp from the ptable and look at whether there is anything between i and ip and jp and j
 * This should ensure that I am always
*/
void generate_pt(std::string &structure, std::vector<int> &fres, std::vector<int> &up, int n){
    std::vector<int> paren,square;
    int count = 0;

    for(int j = 0;j<n;++j){
        switch(structure[j]){
            case '(': paren.push_back(j); count = 0; break;
            case '[': square.push_back(j); count = 0; break;
            case ')': { int i = paren.back();  paren.pop_back();  fres[i] = j; fres[j] = i; count = 0; break; }
            case ']': { int i = square.back(); square.pop_back(); fres[i] = j; fres[j] = i; count = 0; break; }
            default:  up[j] = ++count; break;
        }
    }
    for(auto *stack : {&paren, &square}){
        if(!stack->empty()){
            std::cerr << "Error: unmatched brackets in: " << structure << std::endl;
            exit(1);
        }
    }
}
// i and j are the structured part
static inline bool empty_region(const std::vector<int> &up, int i, int j){
    return up[j-1] >= j-i-1;
}
// By having a function that deals with the fatgraph index on a per bracket type basis, it will be easier to expand this later if needed
static inline void process_bracket(std::string &structure, std::string &fatgraph, std::string &fatgraph_full, cand_pos_t j,const std::vector<int> &fres,const  std::vector<int> &up,std::vector<int> &stack, char open, char close){
    if(structure[j] == open){
        /**
         * Three cases: 1) There is nothing in there yet, so it's new and should be added
         * 2 and 3) There is something (i.e. another bracket type) between the last added bracket of this type
         * and our current one; therefore, it's not and internal loop/stack and should be added
         */
        bool update = stack.empty() || !empty_region(up,stack.back(),j) ||  !empty_region(up,fres[j],fres[stack.back()]);
        if(update){
            fatgraph_full[j] = open;
            fatgraph += open;
        }
        stack.push_back(j);
    } else if(structure[j] == close){ 
        cand_pos_t i = stack.back();
        stack.pop_back();
        // Fatgraph full is the n-length version which makes it easier to do checks
        // If we have it in fatgraph full, then that means it was unique and the close should be added.
        if(fatgraph_full[i] == open){
            fatgraph_full[j] = close;
            fatgraph+=close;
        }
    }
}

std::string generate_fatgraph(std::string &structure,const std::vector<int> &fres,const  std::vector<int> &up, const cand_pos_t n){
    std::vector<int> paren;
    std::vector<int> sb;
    paren.reserve(n/2);
    std::string fatgraph_full = std::string(n,'.');
    std::string fatgraph = "";
   for(cand_pos_t j = 0;j<n;++j){
    process_bracket(structure,fatgraph,fatgraph_full,j,fres,up,paren,'(',')');
    process_bracket(structure,fatgraph,fatgraph_full,j,fres,up,sb,'[',']');
   }
   return fatgraph;
}

/**
 * Write out a consistent format where parentheses come first, then square brackets, then etc.
 */
std::string canonicalize_fatgraph(const std::string &fatgraph){
    std::string canonicalized = fatgraph;
    const std::string openers = "([{<";
    const std::string closers = ")]}>";
    auto opener_idx = [](char c) -> int {
        switch(c){
            case '(': return 0;
            case '[': return 1;
            case '{': return 2;
            case '<': return 3;
            default: return -1;
        }
    };
    std::unordered_map<char,char> remap;
    int next_type = 0;

    for(char c : fatgraph){
        int idx = opener_idx(c);
        if(idx >= 0 && remap.count(c) == 0){
            remap[c]           = openers[next_type];
            remap[closers[idx]] = closers[next_type];
            ++next_type;
        }
    }

    for(char &c : canonicalized){
        if(remap.count(c)) c = remap[c];
    }

    return canonicalized;
}

std::string W_final_pf::get_fatgraph(std::string structure){
    const cand_pos_t n = structure.length();
    std::vector<int> fres(n,-2);
    std::vector<int> up(n,0);
    generate_pt(structure,fres,up,n);

    std::string fatgraph = generate_fatgraph(structure,fres,up,n);
    fatgraph = canonicalize_fatgraph(fatgraph);
    return fatgraph;
}