"~/Nakano_RNAseq/network_analysis/script/ATTEDII/Original_ATTEDPCC_RNASeq.R"
#supple----
T.Exp <- read.table("~/bigdata/yasue/ATTEDII/RNA-Seq/Ath-r.v17-11.G22760-S2120.combat.expression.txt", sep = "\t", quote = "", header = T, row.names = 1, stringsAsFactors = F)
T.ID <- rownames(T.Exp)[order(as.numeric(rownames(T.Exp)))]
#write.table(T.ID,"~/bigdata/yasue/ATTEDII/RNA-Seq/fileID.txt",col.names = F,row.names = F,quote = F)
before <- proc.time()
#puckage----
library(Hmisc)
#input----
T.Description <- read.table("~/bigdata/yasue/ATTEDII/RNA-Seq/Convert_EnterGeneID_to_TAIRID.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#processing data----
T.Exp <- T.Exp[match(T.Description$From, rownames(T.Exp)), ]
rownames(T.Exp) <- T.Description$To[match(rownames(T.Exp), T.Description$From)]
T.Exp <- t(T.Exp)
T.PCC <- rcorr(as.matrix(T.Exp), type = "pearson")
T.PCC <- T.PCC$r
temp <- combn(rownames(T.PCC), 2)
ATTED.PCC <- data.frame(Source = temp[1, ],
                        Target = temp[2, ],
                        PCC = lower.tri(T.PCC),
                        stringsAsFactors = F
                        )
#save
saveRDS(object = ATTED.PCC, file = "~/bigdata/yasue/ATTEDII/RNA-Seq/ATTEDPCC.rds")
#elapsed time----
after <- proc.time()
print(after - before)#1.427 sec
#remove object----
rm(list = ls())