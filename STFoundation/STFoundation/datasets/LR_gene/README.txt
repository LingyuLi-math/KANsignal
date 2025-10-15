## 2025.10.15 Convert LRpair to LRgene
-- Code:
    cd /mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation
    python STFoundation/datasets/LR_gene/LRpair2LRgene.py --species human


## For human
    (KAN) lingyu@GPU2:/mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation$ python STFoundation/datasets/LR_gene/LRpair2LRgene.py --species human
    CellChatDBv2 LRI dim: (3233, 27)
    CellChatDB_v1 Number: 1920
    CellChatDB_v2 Number: 1313
    ['TGFB1' 'TGFB2' 'TGFB3' ... 'PROS1' 'PLAU' 'VCAM1']
    ['TGFbR1_R2' 'TGFbR1_R2' 'TGFbR1_R2' ... 'MERTK' 'PLAUR' 'ITGAD_ITGB2']
    LRgene_df shape: (1404, 1)
    CellChatDB_v1 Gene Number: 963
    CellChatDB_v2 Gene Number: 1404
    Saved LR genes to: LR_gene/LRgene_CellChatDBv2_baseline_human.csv

## For mouse
    (KAN) lingyu@GPU2:/mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation$ python STFoundation/datasets/LR_gene/LRpair2LRgene.py --species mouse
    CellChatDBv2 LRI dim: (3379, 27)
    CellChatDB_v1 Number: 1993
    CellChatDB_v2 Number: 1386
    ['Tgfb1' 'Tgfb2' 'Tgfb3' ... 'Pros1' 'Plau' 'Vcam1']
    ['TGFbR1_R2' 'TGFbR1_R2' 'TGFbR1_R2' ... 'Mertk' 'Plaur' 'ITGAD_ITGB2']
    LRgene_df shape: (1546, 1)
    CellChatDB_v1 Gene Number: 963
    CellChatDB_v2 Gene Number: 1546
    Saved LR genes to: LR_gene/LRgene_CellChatDBv2_baseline_mouse.csv