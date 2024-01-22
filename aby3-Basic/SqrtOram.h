#pragma once
#include "./Oram/include/oram.h"
#include "Shuffle.h"
#include "Basics.h"

// struct boolShare{
//     std::array<bool, 2> bshares;

//     void set_false(){
//         bshares[0] = 1; bshares[1] = 1;
//     }
//     bool reveal(){
//         return bshares[0] ^ bshares[1];
//     }

// };

class ABY3PosMap : PosMap<aby3::si64>{
    public:
        std::vector<boolShare> usage_map;
        ABY3PosMap(size_t len, size_t pack, size_t S, std::vector<aby3::si64>& permutation)
            : PosMap<aby3::si64>(len, pack, S, permutation) {}

        aby3::si64 access(aby3::si64 index, boolShare fake){
            aby3::si64 physical_index; 
            if(linear){
                // bool done = false;
                boolShare done;
                done.set_false();

                for(int i=0; i<this->n; i++){
                    // bool s1 = !fake && (i == index);
                    // bool s2 = fake && (!this->usage_map[i]) && (!done);
                    // if(s1 || s2){
                    //     this->usage_map[i] = true;
                    //     done = true;
                    //     physical_index = this->permutation[i];
                    // }
                    boolShare s1, s2;
                }
            }
            else{
                // first look into the stash.
            //     bool found = false;
            //     size_t h = index / pack;
            //     size_t l = index % pack;
            //     for(int i=0; i<this->t; i++){
            //         if(this->stash[i].logicalIndex == h){
            //             found = true;
            //             physical_index = this->stash[i].packedIndices[l];
            //         }
            //     }
            //     if(!found){
            //         // look into the posMap.
            //         physical_index = this->subPosMap->access(h, fake);
            //     }
            // }
            return physical_index;
        }
    }
};