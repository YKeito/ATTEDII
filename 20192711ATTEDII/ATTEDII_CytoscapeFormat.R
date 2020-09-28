#"~/Nakano_RNAseq/network_analysis/script/ATTEDII/20192711ATTEDII/ATTEDII_CytoscapeFormat.R"
#supple----
path <- "~/bigdata/yasue/ATTEDII/Microarray/Ath-m.v15-08.G20836-S15275.rma.mrgeo.d"
T.file <- list.files(path = "~/bigdata/yasue/ATTEDII/Microarray/Ath-m.v15-08.G20836-S15275.rma.mrgeo.d/")
T.fileID <- T.file[order(as.numeric(T.file))]
#write.table(T.fileID,"~/bigdata/yasue/ATTEDII/Microarray/fileID.txt",col.names = F,row.names = F,quote = F)
before <- proc.time()
#input data----
T.DAVID <- read.table("~/bigdata/yasue/ATTEDII/Microarray/Convert_EnterGeneID_to_TAIRID.txt", sep = "\t", quote = "", header = T, stringsAsFactors = F)
#processing data----
n <- 1
ATTED.PCC <- c()
total <- length(T.fileID)
for (n in 1:total) {
  temp <- read.table(file = paste0(path, "/", T.fileID[n]), header = F, stringsAsFactors = F)
  temp <- temp[order(as.numeric(temp[,1])), ]
  ATTED.PCC <- cbind(ATTED.PCC, as.vector(temp[,3]))
  colnames(ATTED.PCC)[n] <- T.fileID[n]
  print(total-n)
  n <- n+1
}
rownames(ATTED.PCC) <- temp$V1
#saveRDS(object = ATTED.PCC, file = "~/bigdata/yasue/ATTEDII/Microarray/ATTEDPCC.rds")
temp <- combn(rownames(ATTED.PCC), 2)
T.data <- data.frame(Source = temp[1, ],
                     Target = temp[2, ],
                     PCC = ATTED.PCC[lower.tri(ATTED.PCC)],
                     stringsAsFactors = F
                     )
T.data$Source <- T.DAVID$To[match(T.data$Source, T.DAVID$From)]
T.data$Target <- T.DAVID$To[match(T.data$Target, T.DAVID$From)]
T.data <- na.omit(T.data)
#saveRDS(object = T.data, file = "~/bigdata/yasue/ATTEDII/Microarray/ATTEDPCC_CytoscapeFormat.rds")
#elapsed time----
after <- proc.time()
print(after - before)#1.427 sec
#remove object----
rm(list = ls())