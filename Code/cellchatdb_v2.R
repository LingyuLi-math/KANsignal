## 2025.10.15 - 使用python得不到东西，尝试用R看一下


rm(list=ls())
# path = '/home/lingyu/ssd/Python/'
# path = 'E:\\Server\\woody\\Python\\'
# path = 'E:\\Server\\GPU2\\Python\\'

## in my PC
# path = 'D:\\OneDrive - The University Of Hong Kong\\公派出国\\博士后\\HuangYH\\Project\\25.03.17_LRsimulate\\CellChatDB\\CellChat\\data'
# setwd(path)

## in GPU2
path = '/mnt/lingyu/nfs_share2/Python/'
setwd(paste(path, 'STPatch/KANsignal_GPU2/STFoundation/0_CellChatDBv2', sep = ''))
getwd()


## Load CellChat human 
load("CellChatDB.human.rda")  
str(CellChatDB.human, max.level=1)
# List of 4
# $ interaction:'data.frame':	3233 obs. of  28 variables:
#   $ complex    :'data.frame':	338 obs. of  5 variables:
#   $ cofactor   :'data.frame':	32 obs. of  16 variables:
#   $ geneInfo   : tibble [26,827 × 9] (S3: tbl_df\\tbl\\data.frame)


## save

## in my pC
# setwd(paste(path, '\\Processed_human', sep = ''))

## in GPU2
setwd(paste(path, 'STPatch/KANsignal_GPU2/STFoundation/STFoundation/datasets/LR_data/', sep = ''))
getwd()

database = 'human-CellChatv2_'
write.csv(CellChatDB.human$interaction, paste(database, 'interaction.csv', sep = ''))
write.csv(CellChatDB.human$complex, paste(database, 'complex.csv', sep = ''))
write.csv(CellChatDB.human$cofactor, paste(database, 'cofactor.csv', sep = ''))
write.csv(CellChatDB.human$geneInfo, paste(database, 'geneInfo.csv', sep = ''))
