#include "oram.h"

// template <typename T>
// inline void print_vector(std::vector<T>& data){
//     size_t len = data.size();
//     for(size_t i = 0; i < len; i++){
//         std::cout << data[i] << " ";
//     }
//     std::cout << std::endl;
// }

// template <typename T>
// inline void print_vector(std::vector<PackedIndex<T>>& data){
//     size_t len = data.size();
//     for(size_t i = 0; i < len; i++){
//         std::cout << "index: " << data[i].logicalIndex << " ";
//         std::cout << "data: ";
//         for(size_t j = 0; j < data[i].pack; j++){
//             std::cout << data[i].packedIndices[j] << " ";
//         }
//     }
//     std::cout << std::endl;
// }

// template <typename T>
// struct is_packed_index : std::false_type {};

// template <typename T>
// struct is_packed_index<PackedIndex<T>> : std::true_type {};

// template <typename T>
// inline void print_vector_impl(std::vector<T>& data, std::false_type){
//     size_t len = data.size();
//     for(size_t i = 0; i < len; i++){
//         std::cout << data[i] << " ";
//     }
//     std::cout << std::endl;
// }

// template <typename T>
// inline void print_vector_impl(std::vector<PackedIndex<T>>& data, std::true_type){
//     size_t len = data.size();
//     for(size_t i = 0; i < len; i++){
//         std::cout << "index: " << data[i].logicalIndex << " ";
//         std::cout << "data: ";
//         for(size_t j = 0; j < data[i].pack; j++){
//             std::cout << data[i].packedIndices[j] << " ";
//         }
//     }
//     std::cout << std::endl;
// }

// template <typename T>
// inline void print_vector(std::vector<T>& data){
//     print_vector_impl(data, is_packed_index<T>{});
// }

void test_permutation();
void test_pos_map();