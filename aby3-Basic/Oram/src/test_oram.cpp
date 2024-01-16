#include "../include/oram.h"
#include "../include/test.h"

void test_permutation(){
    size_t len = 10;
    std::vector<int> data(len);
    for(int i=0; i<len; i++){
        data[i] = i;
    }
    std::vector<size_t> permutation;
    get_permutation(len, permutation);
    permutate(permutation, data);

    for(int i=0; i<len; i++){
        if(data[permutation[i]] != i){
            std::cout << "permutation error" << std::endl;
            std::cout << "permutation: " << std::endl;
            print_vector(permutation);
            std::cout << "data: " << std::endl;
            print_vector(data);
            exit(1);
        }
    }
}

void test_pos_map(){
    size_t len = 10, pack = 2, S=5;
    std::vector<int> data(len);
    for(int i=0; i<len; i++) data[i] = i;
    std::vector<size_t> permutation;
    get_permutation(len, permutation);
    print_vector(permutation);
    permutate(permutation, data);
    print_vector(data);

    PosMap<size_t> pos_map(len, pack, S, permutation);
    for(int i=0; i<len; i++){
        if(data[pos_map.access(i, false)] != i){
            std::cout << "position map error" << std::endl;
        }
    }
    for(int i=len-1; i>-1; i--){
        if(data[pos_map.access(i, false)] != i){
            std::cout << "position map error" << std::endl;
        }
    }
}

int main(){
    std::srand(16);
    test_permutation();
    test_pos_map();
    return 0;
}