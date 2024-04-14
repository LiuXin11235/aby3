#!/bin/bash

# source_ids=(170 171 172 173) # 139, 185, 24, 60.
source_ids=(171 172 173 174 175 176 177 178) # 139

i=1;
for source_id in ${source_ids[@]}; do
    echo $i | tee /proc/irq/$source_id/smp_affinity_list
    i=$((i+1))
done

declare -A target_irqs
target_irqs[178]="234 179 180 181 182 183 184"
target_irqs[171]="185 186 187 188 189 190 191"
target_irqs[172]="192 193 194 195 196 197 198"
target_irqs[173]="199 200 201 202 203 204 205"
target_irqs[174]="206 207 208 209 210 211 212"    
target_irqs[175]="213 214 215 216 217 218 219"
target_irqs[176]="220 221 222 223 224 225 226"
target_irqs[177]="227 228 229 230 231 232 233"


for source_id in ${source_ids[@]}; do
    read -a target_irqs_list <<< ${target_irqs[$source_id]}
    for target_irq in ${target_irqs_list[@]}; do
        echo "Binding IRQ $target_irq to CPU with the same as $source_id";
        cat /proc/irq/$source_id/smp_affinity_list | tee /proc/irq/$target_irq/smp_affinity_list
    done
done