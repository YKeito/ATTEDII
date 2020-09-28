ATTEDII_pair <- paste0(CoExp_absth04_cyto$souce, CoExp_absth04_cyto$target)
ATTEDII_pair <- data.frame(AGI_pair = mastercluster_pair,
                           PCC = CoExp_absth04_cyto$PCC
                           )

filename <- list.files("~/Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY20/", pattern=".txt", full.names = T)
temp2 <- c("12h", "1h", "24h", "3h")
object_name <- c()
object_name_all <- c()
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("CY20", "_", temp2[n], "_DEGs")
  assign(object_name, read.table(filename[n],  header=T, sep="\t", stringsAsFactors = F))
  n <- n+1
}



AGI_list <- list(rownames(CY15_1h_FDR005), rownames(CY15_3h_FDR005), rownames(CY15_12h_FDR005), rownames(CY15_24h_FDR005),
                  rownames(CY16_1h_FDR005), rownames(CY16_3h_FDR005), rownames(CY16_12h_FDR005), rownames(CY16_24h_FDR005),
                  rownames(CY20_1h_FDR005), rownames(CY20_3h_FDR005), rownames(CY20_12h_FDR005), rownames(CY20_24h_FDR005)
                  )
names(AGI) <- list("CY15_1h_FDR005", "CY15_3h_FDR005", "CY15_12h_FDR005", "CY15_24h_FDR005",
                   "CY16_1h_FDR005", "CY16_3h_FDR005", "CY16_12h_FDR005", "CY16_24h_FDR005",
                   "CY20_1h_FDR005", "CY20_3h_FDR005", "CY20_12h_FDR005", "CY20_24h_FDR005"
                   )

i <- 1
for(i in i:length(AGI_list)){
  CY <- combn(AGI_list[[i]], 2)
  pair <- paste0(CY[1, ], CY[2, ])
  edge <- intersect(mastercluster_pair$AGI_pair, pair)
  PCC <- mastercluster_pair$PCC[match(edge, mastercluster_pair$AGI_pair)]
  souce <- c()
  target <- c()
  m <- 1
  for(m in m:length(edge)){
    souce <- c(souce, substr(edge[m], 1, 9))
    target <- c(target, substr(edge[m], 10, 18))
    m <- m+1
  }
  temp <- data.frame(souce = souce,
                     PCC = PCC,
                     target = target
                     )
  File <- paste0("~/Nakano_RNAseq/network_analysis/base/eachCY_131224h/", names(AGI)[i], ".txt")
  write.table(temp, file = File, append=F, quote = F, sep = "\t", row.names = F)
  print(i)
}


