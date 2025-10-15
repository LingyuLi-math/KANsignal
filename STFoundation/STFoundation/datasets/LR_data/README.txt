## 2025.10.15 CellChatDB v1
(version.1.1.3) (mouse: 2,022 pairs, human: 1,940 pairs)

(base) lingyu@woody:~/ssd/Python/FineST/FineST/FineST/datasets/LR_data$ tree
.
├── human-complex_input_CellChatDB.csv
├── human-interaction_input_CellChatDB.csv.gz
├── mouse-complex_input_CellChatDB.csv
├── mouse-interaction_input_CellChatDB.csv.gz
├── zerafish-complex_input_CellChatDB.csv
└── zerafish-interaction_input_CellChatDB.csv.gz

## 2025.10.14 CellChatDB v2
-- Reference: 
    Jin S, Plikus M V, Nie Q. CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics[J]. Nature protocols, 2025, 20(1): 180-219.
-- Update: 
    dd new literature-supported interactions, including both proteins and nonproteins acting as ligands, 
    leading to a total of ~3,300 interactions for both mouse and human.  
-- Add rda data
    D:\OneDrive - The University Of Hong Kong\公派出国\博士后\HuangYH\Project\25.03.17_LRsimulate\CellChatDB\CellChat\data
-- R code: D:\OneDrive - The University Of Hong Kong\公派出国\博士后\HuangYH\Project\25.03.17_LRsimulate\CellChatDB\CellChat\data\cellchatdb_v2.R

## 2025.10.15 Write Rscript
-- /mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation/STFoundation/datasets/LR_data
-- CellChatDBv2_rda2csv.R
-- Run：
    cd /mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation
    Rscript STFoundation/datasets/LR_data/CellChatDBv2_rda2csv.R human

## For human
    (KAN) lingyu@GPU2:/mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation$ Rscript STFoundation/datasets/LR_data/CellChatDBv2_rda2csv.R human
    Species: human 
    RDA file: CellChatDB.human.rda 
    RDA path: 0_CellChatDBv2/CellChatDB.human.rda 
    Current working directory: /mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation 
    Loaded object: CellChatDB.human 
    List of 4
    $ interaction:'data.frame':    3233 obs. of  28 variables:
    $ complex    :'data.frame':    338 obs. of  5 variables:
    $ cofactor   :'data.frame':    32 obs. of  16 variables:
    $ geneInfo   :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':        26827 obs. of  9 variables:
    Exported CellChatDB tables to STFoundation/datasets/LR_data 

## For Mouse
    (KAN) lingyu@GPU2:/mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation$ Rscript STFoundation/datasets/LR_data/CellChatDBv2_rda2csv.R mouse
    Species: mouse 
    RDA file: CellChatDB.mouse.rda 
    RDA path: 0_CellChatDBv2/CellChatDB.mouse.rda 
    Current working directory: /mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation 
    Loaded object: CellChatDB.mouse 
    List of 4
    $ interaction:'data.frame':    3379 obs. of  28 variables:
    $ complex    :'data.frame':    337 obs. of  5 variables:
    $ cofactor   :'data.frame':    33 obs. of  16 variables:
    $ geneInfo   :'data.frame':    25520 obs. of  9 variables:
    Exported CellChatDB tables to STFoundation/datasets/LR_data 