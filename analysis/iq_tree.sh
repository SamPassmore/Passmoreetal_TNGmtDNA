#!/bin/zsh

## 1. Construct a ML tree
### This command checks the best fit model
### Reconstructs the ML tree using the IQTree search algorithm
### Assess Branch support using ultrfast bootstrap. 
iqtree2 -s ../processed_data/papua_m_daru.fasta -B 1000 -T AUTO > iqtree_papuandaru_Mhaplotype.log

iqtree2 -s ../processed_data/papua_n_daru.fasta -B 1000 -T AUTO > iqtree_papuandaru_Nhaplotype.log

## 2. Apply a partition model
#### This model is used to group genes into 'paritions' that share the same evolutionary characteristics
iqtree2 -s ../processed_data/papua_m_daru.fasta -p ../processed_data/papua_m_daru.fasta.splits.nex -B 1000 -T AUTO