all <- union(CoExp_absth04_cyto$souce, CoExp_absth04_cyto$target)

#######MCLNum#############################################################################################################################################################################
#cluster info
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/eachCY_131224h/ATTEDII/th04/output/CY15", pattern=".csv", full.names = T)
temp2 <- c("12h", "1h", "24h", "3h")
dobject_name <- c()
object_name_all <- c()
MCL_Num <- c()
MCL_Num_all <- c()
BC <- c()
BC_all <- c()
degree <- c()
degree_all <- c()
normalize_degree <- c()
normalize_degree_all <- c()
n <- 1
for(n in 1:length(filename)){
  object_name <- paste0("CY15", "_", temp2[n], "_nodeinfo")
  object_name <- assign(object_name, read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  
  temp1 <- intersect(object_name[, "name"], all)
  temp1 <- match(temp1, all)  
  #MCL_Num <- rep(0, times = length(all))
  #MCL_Num[temp1] <- object_name[, "X__mclCluster"]
  #MCL_Num_all <- cbind(MCL_Num_all, MCL_Num)
  
  BC <- rep(0, times = length(all))
  BC[temp1] <- object_name[, "BetweennessCentrality"]
  BC_all <- cbind(BC_all, BC)
  
  degree <- rep(0, times = length(all))
  normalize_degree <- rep(0, times = length(all))
  degree[temp1] <- object_name[, "Degree"]
  normalize_degree[temp1] <- object_name[, "Degree"]/max(object_name[, "Degree"])
  normalize_degree_all <- cbind(normalize_degree_all, normalize_degree)
  degree_all <- cbind(degree_all, degree)
  n <- n+1
}
#colnames(MCL_Num_all) <- paste0("CY15_", temp2, "_MCLNum")
colnames(BC_all) <- paste0("CY15_", temp2, "_BetweennessCentrality")
colnames(degree_all) <- paste0("CY15_", temp2, "_degree")
colnames(normalize_degree_all) <- paste0("CY15_", temp2, "_normalize_degree")

    
AGI_list <- list(rownames(CY15_1h_FDR005), rownames(CY15_3h_FDR005), rownames(CY15_12h_FDR005), rownames(CY15_24h_FDR005))
DEGs <- c()
DEGs_all <- c()
i <- 1
for(i in i:length(AGI_list)){
  temp <- AGI_list[[i]]
  DEGs <- rep("No", times = length(all))
  temp1 <- intersect(temp, all)
  DEGs[match(temp1, all)] <- "Yes"
  DEGs_all <- cbind(DEGs_all, DEGs)
  i <- i+1
}
colnames(DEGs_all) <- c("CY15_1h_DEGs", "CY15_3h_DEGs", "CY15_12h_DEGs", "CY15_24h_DEGs")

TF <- rep("No", times = length(all))
temp1 <- intersect(TF_family$AGI, all)


TF[match(temp1, all)] <- "Yes"
##########################################################################################
ATTEDII_CY15_Table <- data.frame(AGI = all, 
                                 TF,
                                 DEGs_all,
                                 #MCL_Num_all,
                                 BC_all,
                                 degree_all,
                                 normalize_degree_all,
                                 stringsAsFactors = F
                                 )
  
ATTEDII_CY15_Table <- ATTEDII_CY15_Table[, c(1, 2, 
                                             grep("1h", colnames(ATTEDII_CY15_Table)), 
                                             grep("3h", colnames(ATTEDII_CY15_Table)),
                                             grep("12h", colnames(ATTEDII_CY15_Table)),
                                             grep("24h", colnames(ATTEDII_CY15_Table))
                                             )]
    

write.table(ATTEDII_CY15_Table, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/ATTEDII/ATTEDII_CY15_Table_20180803.txt", append=F, quote = F, sep = "\t", row.names = F)

#1h
ATTEDII_CY15_Table_1hsort <- ATTEDII_CY15_Table[order(ATTEDII_CY15_Table$CY15_1h_normalize_degree, decreasing = T), 
                                                c("AGI", "TF", "CY15_1h_DEGs", colnames(ATTEDII_CY15_Table)[grep("normalize", colnames(ATTEDII_CY15_Table))])]

ATTEDII_CY15_Table_1hDEGs_sort <- ATTEDII_CY15_Table_1hsort[ATTEDII_CY15_Table_1hsort$CY15_1h_DEGs != "No", ]
ATTEDII_CY15_1hTFDEGs_sort <- ATTEDII_CY15_Table_1hDEGs_sort[ATTEDII_CY15_Table_1hDEGs_sort$TF != "No", ]

#3h
ATTEDII_CY15_3hsort <- ATTEDII_CY15_Table[order(ATTEDII_CY15_Table$CY15_3h_normalize_degree, decreasing = T), 
                                         c("AGI", "TF", "CY15_3h_DEGs", colnames(ATTEDII_CY15_Table)[grep("normalize", colnames(ATTEDII_CY15_Table))])]

ATTEDII_CY15_3hDEGs_sort <- ATTEDII_CY15_3hsort[ATTEDII_CY15_3hsort$CY15_3h_DEGs != "No", ]
ATTEDII_CY15_3hTFDEGs_sort <- ATTEDII_CY15_3hDEGs_sort[ATTEDII_CY15_3hDEGs_sort$TF != "No", ]


#12h
ATTEDII_CY15_12hsort <- ATTEDII_CY15_Table[order(ATTEDII_CY15_Table$CY15_12h_normalize_degree, decreasing = T), 
                                          c("AGI", "TF", "CY15_12h_DEGs", colnames(ATTEDII_CY15_Table)[grep("normalize", colnames(ATTEDII_CY15_Table))])]
ATTEDII_CY15_12hDEGs_sort <- ATTEDII_CY15_12hsort[ATTEDII_CY15_12hsort$CY15_12h_DEGs != "No", ]
ATTEDII_CY15_12hTFDEGs_sort <- ATTEDII_CY15_12hDEGs_sort[ATTEDII_CY15_12hDEGs_sort$TF != "No", ]

#24h
ATTEDII_CY15_24hsort <- ATTEDII_CY15_Table[order(ATTEDII_CY15_Table$CY15_24h_normalize_degree, decreasing = T), 
                                          c("AGI", "TF", "CY15_24h_DEGs", colnames(ATTEDII_CY15_Table)[grep("normalize", colnames(ATTEDII_CY15_Table))])]
ATTEDII_CY15_24hDEGs_sort <- ATTEDII_CY15_24hsort[ATTEDII_CY15_24hsort$CY15_24h_DEGs != "No", ]
ATTEDII_CY15_24hTFDEGs_sort <- ATTEDII_CY15_24hDEGs_sort[ATTEDII_CY15_24hDEGs_sort$TF != "No", ]


ATTEDII_CY15_intersection <- ATTEDII_CY15_Table[c(ATTEDII_CY15_Table$CY15_1h_DEGs == "Yes" &
                                                   ATTEDII_CY15_Table$CY15_3h_DEGs == "Yes" &
                                                   ATTEDII_CY15_Table$CY15_12h_DEGs == "Yes" &
                                                   ATTEDII_CY15_Table$CY15_24h_DEGs == "Yes"),
                                               c("AGI", "TF", colnames(ATTEDII_CY15_Table)[grep("normalize", colnames(ATTEDII_CY15_Table))])]

ATTEDII_CY15_intersection1hsort <- ATTEDII_CY15_intersection[order(ATTEDII_CY15_intersection$CY15_1h_normalize_degree, decreasing = T), ]
ATTEDII_CY15_TFDEGsintersection <- ATTEDII_CY15_intersection[ATTEDII_CY15_intersection$TF != "No", ]
ATTEDII_CY15_TFDEGsintersection <- ATTEDII_CY15_TFDEGsintersection[order(ATTEDII_CY15_TFDEGsintersection$CY15_1h_normalize_degree, decreasing = T), ]

write.table(ATTEDII_CY15_1hTFDEGs_sort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/ATTEDII/CY15_TF_DEGs_sort/ATTEDII_CY15_1hTFDEGs_1hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(ATTEDII_CY15_3hTFDEGs_sort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/ATTEDII/CY15_TF_DEGs_sort/ATTEDII_CY15_3hTFDEGs_3hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(ATTEDII_CY15_12hTFDEGs_sort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/ATTEDII/CY15_TF_DEGs_sort/ATTEDII_CY15_12hTFDEGs_12hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(ATTEDII_CY15_24hTFDEGs_sort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/ATTEDII/CY15_TF_DEGs_sort/ATTEDII_CY15_24hTFDEGs_24hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(ATTEDII_CY15_TFDEGsintersection, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/ATTEDII/CY15_TF_DEGs_sort/ATTEDII_CY15_TFDEGs_intersection_1hsort.txt", append=F, quote = F, sep = "\t", row.names = F)