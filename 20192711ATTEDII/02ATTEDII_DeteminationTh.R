"~/UMAP_script/ATTEDII_DeteminationTh.R"

before <- proc.time()
#input data----
T.DAVID <- read.table("~/bigdata/yasue/ATTEDII/Microarray/Convert_EnterGeneID_to_TAIRID.txt", sep = "\t", quote = "", header = T, stringsAsFactors = F)
ATTED.PCC <- readRDS(file = "~/bigdata/yasue/UMAP_Project/RDS/ATTEDPCC_CytoscapeFormat.rds")
#processing data----
ATTED.PCC$Source <- T.DAVID$To[match(ATTED.PCC$Source, T.DAVID$From)]
ATTED.PCC$Target <- T.DAVID$To[match(ATTED.PCC$Target, T.DAVID$From)]
ATTED.PCC <- na.omit(ATTED.PCC)
#DeterminationOfNumEdges
NumEdge <- c()
th <- 0.4
while(th <= 0.5){
  NumEdge <- c(NumEdge, sum(ATTED.PCC$PCC > th))
  th <- th + 0.01
}
print(NumEdge)
ATTED.PCC.th <- ATTED.PCC[ATTED.PCC$PCC > 0.45, ]
#save----
write.table(ATTED.PCC.th, file = "~/bigdata/yasue/UMAP_Project/Table/ATTEDPCC_th045.txt", sep = "\t", quote = F, row.names = F)
saveRDS(object = ATTED.PCC.th, file = "~/bigdata/yasue/UMAP_Project/RDS/ATTEDPCC_th045.rds")
#elapsed time----
after <- proc.time()
print(after - before)#740.14 sec, 12.34 min
#remove object----
rm(list = ls())