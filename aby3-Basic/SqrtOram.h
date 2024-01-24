#pragma once
#include "./Oram/include/oram.h"
#include "Shuffle.h"
#include "Basics.h"

class ABY3PosMap : PosMap<boolIndex>{
    public:
        std::vector<boolShare> usage_map;
        ABY3PosMap(size_t len, size_t pack, size_t S, std::vector<boolIndex>& permutation)
            : PosMap<boolIndex>(len, pack, S, permutation) {
                this->usage_map = std::vector<boolShare>(S);
                for(int i=0; i<S; i++){
                    bool_init_false(0, this->usage_map[i]);
                }
            }

        boolIndex access(boolIndex index, boolShare fake){
            boolIndex physical_index; 
            if(linear){
                
                throw std::runtime_error("Not implemented yet.");
                
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
                throw std::runtime_error("Not implemented yet.");
            }
            return physical_index;
        }
};