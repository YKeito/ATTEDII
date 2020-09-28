#CoExpをまとめるスクリプト(ダウンロードし直したバージョン)
path <- "~/bigdata/yasue/ATTEDII/PCC/negishi/Ath-r.v15-08.G25296-S1401.quantile.mrgeo.d"

#pathのディレクトリに含まれるファイル名をすべてT_filenameに格納する
T_filename <- list.files(path)
T_filename <- T_filename[order(as.numeric(T_filename))]
T_filename_path <- paste0(path, "/", T_filename)

#write.table(T_filename,"~/bigdata/yasue/ATTEDII/PCC/negishi/fileID.txt",col.names = F,row.names = F,quote = F)

#GeneID.txtをDAVID(https://david.ncifcrf.gov/tools.jsp)で読み込んで
#AGIナンバーにコンバートする
T_input <- read.csv("bigdata/yasue/ATTEDII/PCC/negishi/20180802_GeneID.csv", header = T)
T_No <- match(T_filename, T_input$From)
T_AGI <- T_input$To[T_No]

#ATTEDから取得した共発現解析データからPCCの値を抜き出しまとめる
#ファイルに含まれる全GeneIDの組み合わせのPCCをCoExp_PCCに格納する
n <- 1
total <- length(T_filename_path)
CoExp_PCC <- matrix(0, nrow = total, ncol = total, byrow=F)
for (n in 1:total) {
  temp <- read.table(file=T_filename_path[n],header = F)
  temp <- temp[order(temp[,1]),]
  temp <- as.vector(temp[,3])
  CoExp_PCC[, n] <- temp
  print(total - n)
  n <- n+1
}

colnames(CoExp_PCC) <- T_AGI
rownames(CoExp_PCC) <- T_AGI
write.table(CoExp_PCC, "~/bigdata/yasue/ATTEDII/PCC/ATTED_PCC.txt", quote=F, sep="\t", row.names=F, col.names=F)