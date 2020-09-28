#PCC <- CoExp_PCC[lower.tri(CoExp_PCC)]
#temp <- combn(T_AGI, 2)

CoExp_cyto <- data.frame(souce = temp[1, ],
                         PCC = PCC,
                         target = temp[2, ],
                         stringsAsFactors = F)

#write.table(CoExp_cyto, file="~/bigdata/yasue/ATTEDII/PCC/CoExp_cyto_notth.txt", append=F, quote=F, sep="\t", row.names = F, col.names = T)
CoExp_absth06_cyto <- CoExp_cyto[abs(CoExp_cyto$PCC) >  0., ]
write.table(CoExp_absth06_cyto, file="~/bigdata/yasue/ATTEDII/PCC/CoExp_absth06_cyto.txt", append=F, quote=F, sep="\t", row.names = F, col.names = T)